function [grid_out, state_out, R_out, networks_out, ages_out, f_c_out, f_g_out, f_h_out] =...
State_choice(grid_in,c_in,g_in,pH_in,networks_in,ages_in,pro_age_count, q_reduction)
% Function that removes dead cells from grid affter t_n cycles.

N = length(grid_in(:,1));

state = 2*ones(N,N);
R = zeros(N,N);
F = zeros(N,N);
ages = ages_in;

new_grid = grid_in; % temporary grid for calculating cell divisions
shuffleorder = randperm(N*N); % makes sure order is random

f_c = zeros(N,N); % CONSUMPTION/GENERATION RATE OF DIFFERENT NUTRIENTS
f_g = zeros(N,N);
f_h = zeros(N,N);

initial_weights = ANN_loader_met();

    for nr = 1:length(shuffleorder) %     
        %CALCULATES I from remainder,
        %CALCULATES J from floor, 
        i = rem(shuffleorder(nr)-1,N)+1;
        j = floor((shuffleorder(nr)-i)/(N))+1;
            

                    
        if or(grid_in(i,j) == 1, grid_in(i,j) == 2) % if cell is living, in old grid
            
            if i ~= 1 && i ~= N && j ~= 1 && j ~= N % if not on edge
                neighbours = (new_grid(i+1,j) ~= 0)+((new_grid(i-1,j) ~= 0)) + (new_grid(i,j+1) ~= 0) + (new_grid(i,j-1) ~= 0);
            else
                neighbours = 4;
            end

            weights = networks_in{i,j};
            x = [c_in(i,j), g_in(i,j), pH_in(i,j), neighbours]';
            %% KOLLA PÅ R(k+1)
            [state(i,j),R(i,j), MET_TYPE] = ANN_predict(weights,x);
            new_grid(i,j) = MET_TYPE; % SO THAT METABOLISM TYPE CAN CHANGE
            %disp("STATE:");
            %disp(state(i,j,k));
            %F(i,j,k) = mod_consumption( R(i,j,k) );
            %disp(F(i,j,k));
            f_c(i,j) = oxy_rate(    R(i,j), MET_TYPE  );
            f_g(i,j) = glu_rate(    R(i,j), MET_TYPE  );
            f_h(i,j) = hyd_rate(    R(i,j), MET_TYPE  );
            
            if state(i,j) == 0 % if quiescent
                f_c(i,j) = f_c(i,j)/q_reduction;
                f_g(i,j) = f_g(i,j)/q_reduction;
            end
            %disp(MET_TYPE);
 
            killed = 0;
            
            if killed == 0 % NOT KILLED

                if state(i,j) == -1 % if cell chooses apoptosis
                    new_grid(i,j) = 0; % cell dies and disappears
                elseif state(i,j) == 0 % if cell chooses quiescense
                    new_grid(i,j) = MET_TYPE; % cell remains
                elseif state(i,j) == 1 && neighbours < 4 && ages(i,j) >= pro_age_count % if cell proliferates
                    new_grid(i,j) = MET_TYPE; % cell remains;
                    networks_in{i,j} = ANN_mutate(weights);

                    empty_grids = []; %Checks empty spots
                    if new_grid(i,j+1) == 0 % right
                        empty_grids(end+1) = 1;
                    end
                    if new_grid(i+1,j) == 0 % up
                        empty_grids(end+1) = 2;
                    end
                    if new_grid(i,j-1) == 0 % left
                        empty_grids(end+1) = 3;
                    end
                    if new_grid(i-1,j) == 0 % down
                        empty_grids(end+1) = 4;
                    end
                    
                    
                    rand_grid = empty_grids(randsample(length(empty_grids),1));  
                    if rand_grid == 1 % right
                        new_grid(i,j+1) = MET_TYPE; % copy of itself
                        networks_in{i,j+1} = networks_in{i,j};
                        ages(i,j) = 0;
                        ages(i,j+1) = 0;
                    elseif rand_grid == 2 % up
                        new_grid(i+1,j) = MET_TYPE; % copy of itself
                        networks_in{i+1,j} = networks_in{i,j};
                        ages(i,j) = 0;
                        ages(i+1,j) = 0;
                    elseif rand_grid == 3 % left
                        new_grid(i,j-1) = MET_TYPE; % copy of itself
                        networks_in{i,j-1} = networks_in{i,j};
                        ages(i,j) = 0;
                        ages(i,j-1) = 0;
                    elseif rand_grid == 4 % down
                        new_grid(i-1,j) = MET_TYPE; % copy of itself
                        networks_in{i-1,j} = networks_in{i,j};
                        ages(i,j) = 0;
                        ages(i-1,j) = 0;
                    end
                        
                    % NEighbours != 4 here

                else % happens if state == 1 and neighbours = 4
                    % SAME BEHAVIOR AS QUIESCENT
                    new_grid(i,j) = MET_TYPE;  
                end % END IF STATE == 1,0,-1       
            end % END IF NOT KILLED
        end  % END IF NOT ON EDGE
        
        if grid_in(i,j) == 3 % if healthy cell
            networks_in{i,j} = initial_weights; % initial weigths for helathy cells (non-mutated)
            
            if i ~= 1 && i ~= N && j ~= 1 && j ~= N % if not on edge
                neighbours = (new_grid(i+1,j) ~= 0)+((new_grid(i-1,j) ~= 0)) + (new_grid(i,j+1) ~= 0) + (new_grid(i,j-1) ~= 0);
            else
                neighbours = 4;
            end
            
            x = [c_in(i,j), g_in(i,j), pH_in(i,j), neighbours]';
            [state(i,j),R(i,j),MET_TYPE] = ANN_predict(initial_weights,x);
            
            if state(i,j) == -1
               new_grid(i,j) = 0; % cell dies
            else
                state(i,j) = 3; % denotes living healthy cell
            end
            
            f_c(i,j) = oxy_rate(    R(i,j), MET_TYPE  )/q_reduction;
            f_g(i,j) = glu_rate(    R(i,j), MET_TYPE  )/q_reduction;
            f_h(i,j) = hyd_rate(    R(i,j), MET_TYPE  )/q_reduction;
            
        end
        
    end %END FOR NR LOOP
    
 grid_out = new_grid;
 state_out = state;
 R_out = R;
 networks_out = networks_in;
 ages_out = ages;
 f_c_out = f_c;
 f_g_out = f_g;
 f_h_out = f_h;
end