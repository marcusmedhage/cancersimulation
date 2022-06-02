%%% Main code for project
% Requires Diff_oxygen.m, Diff_glucose.m and Diff_hydrogen.m, ANN_loader.m,
% ANN_predict.m, sigmoid.m and inputProcess.m
clear; close all; clc;

tic;

%% Random seed

rng(29);

%% Grid size, cycles

N = 6; % number of cells
T = 1000; % number of cycles
num_of_plots = 1; % How often the cell cycles are plotted
pro_age_count = 10; % iterations before proliferation
t_step = 0.1; % time step for network response
t_n = 5*10; % cycles before a dead cell is remove
q_reduction = 2; % how much metabolism is reduced when quiescent

%% Empty matrices to be filled

c_disc = zeros(N+2,N+2,T+2);
g_disc = zeros(N+2,N+2,T+2);
h_disc = zeros(N+2,N+2,T+2);
pH_disc = zeros(N+2, N+2, T+2);

grid1 = zeros(N+2,N+2,T+2); % Grids where different cells are placed. 
grid2 = zeros(N+2,N+2,T+2);
grid_tot = zeros(N+2,N+2,T+2);

%% Constants

tau = 16*60*60; % cell cycle time in seconds(16 h)
h_size = 25*10^-4; % room step size in cm (real cell size)
n_0 = h_size^(-2); % cell density
c_0 = 1.7*10^-8; % background concentration oxygen
g_0 = 1.3*10^-8; % background concentration glucose
h_0 = 1.0*10^-13; % background concentration hydrogen
%h_0 = 10^(-7.4); % background conc. hydrogen calculated from pH

% Apoptosis threshholds
c_ap = 2.5*10^-9/c_0; % apoptosis threshhold for oxygen
g_ap = 6.5*10^-9/g_0; % apoptosis threshhold for glucose
pH_ap = 6.5; % apoptosis threshhold for pH-level
%h_ap = 10^(-6.5)/(h_0/n_0); % apoptosis threshhold for oxygen
c_m = 0; % anaerobic threshold. Test value for now

%consumption rates
r_c = 2.3*10^-16*tau*n_0/c_0; % base oxygen consumption
r_g_a = 3.8*10^-17*tau*n_0/g_0; % base glucose consuption (aerobic)
r_g_an = 6.9*10^-16*tau*n_0/g_0; % base glucose consuption (anaerobic)
r_h = 1.5*10^-18*tau*n_0/h_0; % base hydrogen increase per cycle


%% Initial conditions

c_in = 0.6*ones(N+2,N+2); % initial conditions oxygen
g_in = 1*ones(N+2,N+2); % initial conditions glucose
h_in = ones(N+2,N+2); % initial conditions hydrogen
%pH_in = -log10(h_in*h_0); % initial conditions pH
pH_in = -log10(h_0*(h_in*10*10/(h_size*0.1)));

c_disc(:,:,1) = c_in;
g_disc(:,:,1) = g_in;
h_disc(:,:,1) = h_in;
pH_disc(:,:,1) = pH_in;

f_c1 = zeros(N+2,N+2); f_c2 = zeros(N+2,N+2);
f_g1 = zeros(N+2,N+2); f_g2 = zeros(N+2,N+2);
f_h1 = zeros(N+2,N+2); f_h2 = zeros(N+2,N+2);

% for grid: -2 = edge -1 = dead cell, 0 = empty cell, 1 = aerobic cell, 2 = anaerobic cell, 3 = healthy cell

grid1(:,:,1) = 3; % healty cells everywhere
% NO GRID2 update

%grid1(:,:,1) = ones(N+2,N+2); %inital conditions (full grid as test)
%grid1(:,:,1) = randi([0 1], N+2,N+2); %inital conditions (random grid as test)
%grid1((N+2)/2:(N+2)/2 + 1, (N+2)/2:(N+2)/2 + 1, 1) = 1; % intial condition (4 cells in middle)
grid1((N+2)/2-1:(N+2)/2 + 2, (N+2)/2-1:(N+2)/2 + 2, 1) = 1; % intial condition (16 cells in middle)
%grid1((N+2)/2,(N+2)/2,1) = 1; % initial conditions (one cell in middle)

grid1(1,:,:) = -2; grid1(N+2,:,:) = -2; grid1(:,1,:) = -2; grid1(:,N+2,:) = -2; % edge
grid2(1,:,:) = -2; grid2(N+2,:,:) = -2; grid2(:,1,:) = -2; grid2(:,N+2,:) = -2; % edge


%% Neural networks

% Vectors to be filled
R1 = zeros(N+2,N+2,T+1); R2 = zeros(N+2,N+2,T+1);
F1 = zeros(N+2,N+2,T+1); F2 = zeros(N+2,N+2,T+1);
state1 = 2*ones(N+2,N+2,T+1);
state1(1,:,:) = -2; state1(N+2,:,:) = -2; state1(:,1,:) = -2; state1(:,N+2,:) = -2; % edge
state2 = 2*ones(N+2,N+2,T+1);
state2(1,:,:) = -2; state2(N+2,:,:) = -2; state2(:,1,:) = -2; state2(:,N+2,:) = -2; % edge

rate_AnAe = zeros(T,1);
total_cancer_cells = zeros(T,1);

% for state: -2 = edge -1 = apoptosis, 0 = quiescence, 1 = proliferation, 2 = empty cell

% Initial networks
ages1 = zeros(N+2,N+2);
ages2 = zeros(N+2,N+2);

weights = ANN_loader_pressure(); % weights for neural network matrices, with metabolism node.
initial_weights = ANN_loader_pressure();

networks1 = cell((N+2),(N+2));
networks2 = cell(N+2,N+2);
for i =1:N+2 % Copies network from initial grid
    for j =1:N+2
        if or(grid1(i,j,1) == 1, grid1(i,j,1) == 2) % if cell is alive
            networks1{i,j} = weights; % intial conditions
            ages1(i,j) = floor(pro_age_count*rand());
        end
    end
end

%% Time stepping

% WE NEED TO PUT CELL CONSUMPTION OUTSIDE HERE.
wait_timer = waitbar(0/T, "RUNNING SIMULATIONS");


for k = 1:T+1
   
    ages1 = ages1 +1;
    ages2 = ages2 +1;
    
    waitbar(k/T, wait_timer);

    %grid(:,:,k) = Dead_cell_removal(grid(:,:,1:k) ,t_n, k); % removes dead cells after t_n cycles
    
%   grid1(:,:,k) = Hard_limits_check(grid1(:,:,k),c_disc(:,:,k),g_disc(:,:,k),pH_disc(:,:,k)); % kills cells that exceed hard limit
%  grid2(:,:,k) = Hard_limits_check(grid2(:,:,k),c_disc(:,:,k),g_disc(:,:,k),pH_disc(:,:,k)); % kills cells that exceed hard limit

    % NOTE, CHECK THAT TIME INDEX IS RIGHT, SEEMS WEIRD that k+1 and k.

    
    new_grid1 = grid1(:,:,k); % temporary grid for calculating cell divisions
    new_grid2 = grid2(:,:,k);
    shuffleorder = randperm(2*(N+2)*(N+2)); % makes sure order is random

    f_c = zeros(N+2,N+2); % CONSUMPTION/GENERATION RATE OF DIFFERENT NUTRIENTS
    f_g = zeros(N+2,N+2);
    f_h = zeros(N+2,N+2);

    f_c1 = zeros(N+2,N+2);  f_c2 = zeros(N+2,N+2);
    f_g1 = zeros(N+2,N+2);  f_g2 = zeros(N+2,N+2);
    f_h1 = zeros(N+2,N+2);  f_h2 = zeros(N+2,N+2);

    for nr = 1:length(shuffleorder) %     
        %CALCULATES I from remainder,
        %CALCULATES J from floor, 
        i = rem(shuffleorder(nr)-1,N+2)+1;
        j = floor((shuffleorder(nr)-i)/(N+2)) + 1;
        if j > N+2
            j = j - (N+2);
        end
        layer = floor(shuffleorder(nr) / (length(shuffleorder)/2) ) +1;
%% GRID 1
        if layer == 1 % GRID 1 activated

            if or(grid1(i,j,k) == 1, grid1(i,j,k) == 2) % if cell is living, in old grid
                if i ~= 1 && i ~= N+2 && j ~= 1 && j ~= N+2 % if not on edge
                    neighbours = (new_grid2(i+1,j) ~= 0)+((new_grid2(i-1,j) ~= 0)) + (new_grid2(i,j+1) ~= 0)...
                                + (new_grid2(i,j-1) ~= 0) + (new_grid1(i+1,j) ~= 0) + ...
                                ((new_grid1(i-1,j) ~= 0)) + (new_grid1(i,j+1) ~= 0) + (new_grid1(i,j-1) ~= 0) + (new_grid2(i,j) ~= 0);
                else
                    neighbours = 9;
                end

                weights = networks1{i,j};
                x = [c_in(i,j), g_in(i,j), pH_in(i,j), neighbours]';
                %% KOLLA PÃ… R(k+1)
                [state1(i,j,k),R1(i,j,k), MET_TYPE] = ANN_predict(weights,x);
%                state1(i,j,k) = 1; % FORCES PROLIFERATION
                new_grid1(i,j) = MET_TYPE; % SO THAT METABOLISM TYPE CAN CHANGE
                f_c1(i,j) = oxy_rate(    R1(i,j,k), MET_TYPE  );
                f_g1(i,j) = glu_rate(    R1(i,j,k), MET_TYPE  );
                f_h1(i,j) = hyd_rate(    R1(i,j,k), MET_TYPE  );
                if state1(i,j,k) == 0 % if quiescent
                    f_c1(i,j) = f_c1(i,j)/q_reduction;
                    f_g1(i,j) = f_g1(i,j)/q_reduction;
                    f_h1(i,j) = f_h1(i,j)/q_reduction;
                    ages1(i,j) = 0;
                end

            if state1(i,j,k) == -1 % if cell chooses apoptosis
                new_grid1(i,j) = 0; % cell dies and disappears
            elseif state1(i,j,k) == 0 % if cell chooses quiescense
                new_grid1(i,j) = MET_TYPE; % cell remains
%% PROLIFERATION OCCURS HERE, need to rewrite this, retrain, network, it gets strange when we have twice the number of neighbours
            elseif state1(i,j,k) == 1 && neighbours < 9 && ages1(i,j) >= pro_age_count % if cell proliferates
                new_grid1(i,j) = MET_TYPE; % cell remains;

                empty_grids = []; %Checks empty spots, and if pressure sufficient,
                % CHECK EXACTLY SAME POINT ON BOTH
                if (new_grid1(i,j+1) == 0 || new_grid2(i,j+1) == 0) && new_grid2(i,j+1) ~= 3 % right
                    empty_grids(end+1) = 1;
                end
                if (new_grid1(i+1,j) == 0 || new_grid2(i+1,j) == 0) && new_grid2(i+1,j) ~= 3% up
                    empty_grids(end+1) = 2;
                end
                if (new_grid1(i,j-1) == 0 || new_grid2(i,j-1) == 0) && new_grid2(i,j-1) ~= 3% left
                    empty_grids(end+1) = 3;
                end
                if (new_grid1(i-1,j) == 0 || new_grid2(i-1,j) == 0) && new_grid2(i-1,j) ~= 3% down
                    empty_grids(end+1) = 4;
                end 
                if new_grid1(i,j) == 0 || new_grid2(i,j) == 0  && new_grid2(i,j) ~= 3 %NOTE, order of grid 1 or 2 irrelevant, tries to place in grid 1 first
                    empty_grids(end+1) = 5; % NO PRESSURE CHECK IF SAME.
                end % The case if there is room in same square as before.
                % MÅSTE UPPDATERAS MED CHECK FÖR VILKA nya grids som är
                % lediga, samt vilka tryck som finns och om
                % tillräckligt stort 
                % HERE IT IS NEEDED TO PLACE IN EITHER GRID 1 or two,
                % possible to try grid 1 first and check and if
                % occupied place in grid 2.
                if isempty(empty_grids) == 1 % if prohibited to grow due to healthy cell in grid 2
                        rand_grid = -1;
                else
                    rand_grid = empty_grids(randsample(length(empty_grids),1)); 
                    network = ANN_mutate(weights);
                    networks1{i,j} = network;
                end


                if rand_grid == 1 % right
                    if new_grid1(i,j+1) == 0
                        new_grid1(i,j+1) = MET_TYPE; % copy of itself
                        networks1{i,j+1} = network;
                        ages1(i,j) = 0;
                        ages1(i,j+1) = 0;
                    elseif new_grid2(i,j+1) == 0 && new_grid1(i,j+1) ~= 3 % ELSE placed in newgrid 2
                        new_grid2(i,j+1) = MET_TYPE; % copy of itself
                        networks2{i,j+1} = network;
                        ages2(i,j) = 0;
                        ages2(i,j+1) = 0;
                    end

                elseif rand_grid == 2 % up
                    if new_grid1(i+1,j) == 0
                        new_grid1(i+1,j) = MET_TYPE; % copy of itself
                        networks1{i+1,j} = network;
                        ages1(i,j) = 0;
                        ages1(i+1,j) = 0;
                    elseif new_grid2(i+1,j) == 0 && new_grid1(i+1,j) ~= 3
                        new_grid2(i+1,j) = MET_TYPE; % copy of itself
                        networks2{i+1,j} = network;
                        ages2(i,j) = 0;
                        ages2(i+1,j) = 0;
                    end
                elseif rand_grid == 3 % left
                    if new_grid1(i,j-1) == 0
                        new_grid1(i,j-1) = MET_TYPE; % copy of itself
                        networks1{i,j-1} = networks1{i,j};
                        ages1(i,j) = 0;
                        ages1(i,j-1) = 0;
                    elseif new_grid2(i,j-1) == 0 && new_grid1(i,j-1) ~= 3
                        new_grid2(i,j-1) = MET_TYPE; % copy of itself
                        networks2{i,j-1} = networks1{i,j};
                        ages2(i,j) = 0;
                        ages2(i,j-1) = 0;
                    end
                elseif rand_grid == 4 % down
                    if new_grid1(i-1,j) == 0
                        new_grid1(i-1,j) = MET_TYPE; % copy of itself
                        networks1{i-1,j} = networks1{i,j};
                        ages1(i,j) = 0;
                        ages1(i-1,j) = 0;
                    elseif new_grid2(i-1,j) == 0 && new_grid1(i-1,j) ~= 3
                        new_grid2(i-1,j) = MET_TYPE; % copy of itself
                        networks2{i-1,j} = networks1{i,j};
                        ages2(i,j) = 0;
                        ages2(i-1,j) = 0;
                    end
                elseif rand_grid == 5 % RANDGRID == 5, 
                   if new_grid1(i,j) == 0 %IF is in 2 and places in 1
                        new_grid1(i,j) = MET_TYPE; % copy of itself
                        networks1{i,j} = networks2{i,j};
                        ages1(i,j) = 0;
                        ages2(i,j) = 0;
                    elseif new_grid2(i,j) == 0 && new_grid1(i,j) ~= 3
                        new_grid2(i,j) = MET_TYPE; % copy of itself
                        networks2{i,j} = networks1{i,j};
                        ages1(i,j) = 0;
                        ages2(i,j) = 0;
                    end 

                end % END RANDGRID CHECK

                % NEighbours != 4 here

            else % happens if state == 1 and neighbours = 9
                % SAME BEHAVIOR AS QUIESCENT
                new_grid1(i,j) = MET_TYPE;  
            end % END IF STATE == 1,0,-1       
        end  % END IF NOT ON EDGE
        
        if grid1(i,j,k) == 3 % if healthy cell
            
            if i ~= 1 && i ~= N+2 && j ~= 1 && j ~= N+2 % if not on edge
                neighbours = (new_grid1(i+1,j) ~= 0)+((new_grid1(i-1,j) ~= 0)) + (new_grid1(i,j+1) ~= 0) + (new_grid1(i,j-1) ~= 0);
            else
                neighbours = 4;
            end
            
            x = [c_in(i,j), g_in(i,j), pH_in(i,j), neighbours]';
            state1(i,j,k) = cell_resp(c_in(i,j), g_in(i,j), pH_in(i,j), neighbours);
            
            if state1(i,j,k) == -1
               new_grid1(i,j) = 0; % cell dies
            else
                state1(i,j,k) = 3; % denotes living healthy cell
            end
            
            
%             f_c1(i,j) = oxy_rate(    R1(i,j), MET_TYPE  )/q_reduction;
%             f_g1(i,j) = glu_rate(    R1(i,j), MET_TYPE  )/q_reduction;
%             f_h1(i,j) = hyd_rate(    R1(i,j), MET_TYPE  )/q_reduction;
            % Healthy cells consume nothing (equlibrium)
            f_c1(i,j) = 0;
            f_g1(i,j) = 0;
            f_h1(i,j) = 0;
            
        end
%% LAYER 2
        elseif layer == 2 % ELSE LAYER 2  
            if or(grid2(i,j,k) == 1, grid2(i,j,k) == 2) % if cell is living, in old grid
                    if i ~= 1 && i ~= N+2 && j ~= 1 && j ~= N+2 % if not on edge
                        neighbours = (new_grid2(i+1,j) ~= 0)+((new_grid2(i-1,j) ~= 0)) + ...
                       (new_grid2(i,j+1) ~= 0) + (new_grid2(i,j-1) ~= 0) + (new_grid1(i+1,j) ~= 0) ...
                       +((new_grid1(i-1,j) ~= 0)) + (new_grid1(i,j+1) ~= 0) + (new_grid1(i,j-1) ~= 0) + (new_grid1(i,j) ~= 0);
                    else 
                        neighbours = 9;
                    end
                    %% CHANGE THIS AND RETRAIN NETWORK
                    weights = networks2{i,j};
                    x = [c_in(i,j), g_in(i,j), pH_in(i,j), neighbours]';
                    %% KOLLA PÃ… R(k+1)
                    [state2(i,j,k),R2(i,j,k), MET_TYPE] = ANN_predict(weights,x);
%                    state2(i,j,k) = 1; % FORCES PROLIFERATION
                    new_grid2(i,j) = MET_TYPE; % SO THAT METABOLISM TYPE CAN CHANGE
                    f_c2(i,j) = oxy_rate(    R2(i,j,k), MET_TYPE  );
                    f_g2(i,j) = glu_rate(    R2(i,j,k), MET_TYPE  );
                    f_h2(i,j) = hyd_rate(    R2(i,j,k), MET_TYPE  );
                    if state2(i,j) == 0 % if quiescent
                        f_c2(i,j) = f_c2(i,j)/q_reduction;
                        f_g2(i,j) = f_g2(i,j)/q_reduction;
                        f_h2(i,j) = f_h2(i,j)/q_reduction;
                        ages2(i,j) = 0;
                    end

                if state2(i,j,k) == -1 % if cell chooses apoptosis
                    new_grid2(i,j) = 0; % cell dies and disappears
                elseif state2(i,j,k) == 0 % if cell chooses quiescense
                    new_grid2(i,j) = MET_TYPE; % cell remains
    %% PROLIFERATION OCCURS HERE, need to rewrite this, retrain, network, it gets strange when we have twice the number of neighbours
                elseif state2(i,j,k) == 1 && neighbours < 9 && ages2(i,j) >= pro_age_count % if cell proliferates
                    new_grid2(i,j) = MET_TYPE; % cell remains;
                    network = ANN_mutate(weights);
                    networks2{i,j} = network;

                    empty_grids = []; %Checks empty spots, and if pressure sufficient,
                    % CHECK EXACTLY SAME POINT ON BOTH
                    if (new_grid1(i,j+1) == 0 || new_grid2(i,j+1) == 0) && new_grid1(i,j+1) ~= 3 % right
                        empty_grids(end+1) = 1;
                    end
                    if (new_grid1(i+1,j) == 0 || new_grid2(i+1,j) == 0) && new_grid1(i+1,j) ~= 3% up
                        empty_grids(end+1) = 2;
                    end
                    if (new_grid1(i,j-1) == 0 || new_grid2(i,j-1) == 0) && new_grid1(i,j-1) ~= 3% left
                        empty_grids(end+1) = 3;
                    end
                    if (new_grid1(i-1,j) == 0 || new_grid2(i-1,j) == 0) && new_grid1(i-1,j) ~= 3% down
                        empty_grids(end+1) = 4;
                    end 
                    if new_grid1(i,j) == 0 || new_grid2(i,j) == 0  && new_grid1(i,j) ~= 3 %NOTE, order of grid 1 or 2 irrelevant, tries to place in grid 1 first
                        empty_grids(end+1) = 5; % NO PRESSURE CHECK IF SAME.
                    end % The case if there is room in same square as before.
                    % MÅSTE UPPDATERAS MED CHECK FÖR VILKA nya grids som är
                    % lediga, samt vilka tryck som finns och om
                    % tillräckligt stort
                    if isempty(empty_grids) == 1 % if prohibited to grow due to healthy cell in grid1
                        rand_grid = -1;
                    else
                       rand_grid = empty_grids(randsample(length(empty_grids),1)); 
                       network = ANN_mutate(weights);
                       networks2{i,j} = network;
                    end 
                    % HERE IT IS NEEDED TO PLACE IN EITHER GRID 1 or two,
                    % possible to try grid 1 first and check and if
                    % occupied place in grid 2. 
                    
                    if rand_grid == 1 % right
                        if new_grid1(i,j+1) == 0 && new_grid2(i,j+1) ~= 3
                            new_grid1(i,j+1) = MET_TYPE; % copy of itself
                            networks1{i,j+1} = network;
                            ages1(i,j) = 0;
                            ages1(i,j+1) = 0;
                        else % ELSE placed in newgrid 2
                            new_grid2(i,j+1) = MET_TYPE; % copy of itself
                            networks2{i,j+1} = network;
                            ages2(i,j) = 0;
                            ages2(i,j+1) = 0;
                        end

                    elseif rand_grid == 2 % up
                        if new_grid1(i+1,j) == 0  && new_grid2(i+1,j) ~= 3
                            new_grid1(i+1,j) = MET_TYPE; % copy of itself
                            networks1{i+1,j} = network;
                            ages1(i,j) = 0;
                            ages1(i+1,j) = 0;
                        else
                            new_grid2(i+1,j) = MET_TYPE; % copy of itself
                            networks2{i+1,j} = network;
                            ages2(i,j) = 0;
                            ages2(i+1,j) = 0;
                        end
                    elseif rand_grid == 3 % left
                        if new_grid1(i,j-1) == 0  && new_grid2(i,j-1) ~= 3
                            new_grid1(i,j-1) = MET_TYPE; % copy of itself
                            networks1{i,j-1} = network;
                            ages1(i,j) = 0;
                            ages1(i,j-1) = 0;
                        else % layer 2
                            new_grid2(i,j-1) = MET_TYPE; % copy of itself
                            networks2{i,j-1} = network;
                            ages2(i,j) = 0;
                            ages2(i,j-1) = 0;
                        end
                    elseif rand_grid == 4 % down
                        if new_grid1(i-1,j) == 0  && new_grid2(i-1,j) ~= 3
                            new_grid1(i-1,j) = MET_TYPE; % copy of itself
                            networks1{i-1,j} = network;
                            ages1(i,j) = 0;
                            ages1(i-1,j) = 0;
                        else
                            new_grid2(i-1,j) = MET_TYPE; % copy of itself
                            networks2{i-1,j} = network;
                            ages2(i,j) = 0;
                            ages2(i-1,j) = 0;
                        end
                    elseif rand_grid == 5 % RANDGRID == 5, 
                       if new_grid1(i,j) == 0  && new_grid2(i,j) ~= 3%IF is in 2 and places in 1
                            new_grid1(i,j) = MET_TYPE; % copy of itself
                            networks1{i,j} = network;
                            ages1(i,j) = 0;
                            ages2(i,j) = 0;
                        else % if is in 1 and places in 2
                            new_grid2(i,j) = MET_TYPE; % copy of itself
                            networks2{i,j} = network;
                            ages1(i,j) = 0;
                            ages2(i,j) = 0;
                        end 

                    end % END RANDGRID CHECK


                    % NEighbours != 4 here

                else % happens if state == 1 and neighbours = 9
                    % SAME BEHAVIOR AS QUIESCENT
                    new_grid2(i,j) = MET_TYPE;  
                end % END IF STATE == 1,0,-1       
            end  % END IF NOT ON EDGE
            

            if grid2(i,j,k) == 3 % if healthy cell
                if i ~= 1 && i ~= N+2 && j ~= 1 && j ~= N+2 % if not on edge
                    neighbours = (new_grid2(i+1,j) ~= 0)+((new_grid2(i-1,j) ~= 0)) + (new_grid2(i,j+1) ~= 0) + (new_grid2(i,j-1) ~= 0);
                else
                    neighbours = 4;
                end

                x = [c_in(i,j), g_in(i,j), pH_in(i,j), neighbours]';
                state2(i,j,k) = cell_resp(c_in(i,j), g_in(i,j), pH_in(i,j), neighbours);

                if state2(i,j,k) == -1
                   new_grid2(i,j) = 0; % cell dies
                else
                    state2(i,j,k) = 3; % denotes living healthy cell
                end

%                 f_c2(i,j) = oxy_rate(    R2(i,j,k), MET_TYPE  )/q_reduction;
%                 f_g2(i,j) = glu_rate(    R2(i,j,k), MET_TYPE  )/q_reduction;
%                 f_h2(i,j) = hyd_rate(    R2(i,j,k), MET_TYPE  )/q_reduction;
                
                % Healthy cells consume nothing (equlibrium)
                f_c1(i,j) = 0;
                f_g1(i,j) = 0;
                f_h1(i,j) = 0;

            end % END HEALTHY CELLS


        end % END IF LAYER 1 ELSE LAYER 2
    end % END SHUFFLEORDER NR LOOP

    f_c = f_c1 + f_c2;
    f_g = f_g1 + f_g2;
    f_h = f_h1 + f_h2;

    grid1(:,:,k+1) = new_grid1;
    grid2(:,:,k+1) = new_grid2;
    
    grid1(:,:,k+1) = Hard_limits_check(grid1(:,:,k+1),c_disc(:,:,k),g_disc(:,:,k),pH_disc(:,:,k)); % kills cells that exceed hard limit
    grid2(:,:,k+1) = Hard_limits_check(grid2(:,:,k+1),c_disc(:,:,k),g_disc(:,:,k),pH_disc(:,:,k)); % kills cells that exceed hard limit
    
    An_count = sum(grid1(:,:,k) == 2) + sum(grid2(:,:,k) == 2);
    Ae_count = sum(grid1(:,:,k) == 1) + sum(grid1(:,:,k) == 1);
    Nec_count = sum(grid1(:,:,k) == -1) + sum(grid1(:,:,k) == -1);
    rate_AnAe(k) = An_count /(Ae_count + An_count);
    total_cancer_cells(k) = sum(Ae_count + An_count + Nec_count);

    disp("A");
    c_in = Diff_oxygen(c_in, f_c , t_step);
    c_disc(:,:,k+1) = c_in;
    
    g_in = Diff_glucose(g_in, f_g, t_step);
    g_disc(:,:,k+1) = g_in;
    
    h_in = Diff_hydrogen(h_in, f_h, t_step);
    h_disc(:,:,k+1) = h_in;
    
    %pH_in = -log10(h_0*h_in);
    pH_in = -log10(h_0*(h_in*10*10/(h_size*0.1)));
    pH_disc(:,:,k+1)  = pH_in;
    
    grid_tot(:,:,k) = grid_total(grid1(:,:,k),grid2(:,:,k));

end

closewaittimer = findall(0,'type', 'figure', "tag",'TMWWaitbar');
delete(closewaittimer);

fprintf ( 'Computation took % g sec \n ' , toc );

%% Visualisation

tic;

visualize_spread(c_disc, g_disc, h_disc, pH_disc, grid1, grid2, grid_tot, state1, state2, num_of_plots,t_step); % plots ocygen, glucose, hydrogen, pH and grid, state

figure(10);
plot(rate_AnAe);
title("Rate of anaerobic cells versus total cells");
xlabel("Timestep");
ylabel("Rate of anaerobic cells");

disp("RATE AN AE: " + rate_AnAe(end));

% figure(11);
% plot(total_cancer_cells);
% title("Number of total cells in tumour");
% xlabel("Timestep");
% ylabel("number of cells cells");


fprintf ( 'Visualisation took % g sec \n ' , toc );