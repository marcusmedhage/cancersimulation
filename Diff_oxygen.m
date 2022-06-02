function c = Diff_oxygen(c_in, f_c, t_end)
% Function that control the oxygen concentration after 1 cell cycle, given
% input oxygen matrix c_in, number of cells N, input grid and conumption R.
h_size = 25*10^-4; % room step size in cm (real cell size)
tau = 16*60*60; % cell cycle time in seconds(16 h)
t_size = 0.001/tau; % time step size in seconds % cell cycles, / (s / cycles)
T = t_end; % end time (1 cycle)
t = 0:t_size:T; % time discretization in cell time base
Ts = length(t);

L = 1;
n_0 = h_size^(-2);
c_0 = 1.7*10^-8;
N = size(c_in, 1)-2;


D_c = 1.8*10^-5*tau/L^2; % Oxygen diffusion constant
r_c = 2.3*10^-16*tau*n_0/c_0; % base oxygen consumption
k_r = 6;
T_r = 0.675;

%Euler forward
c = c_in; % intitial conditions

for k = 2:Ts
    for j = 2:N+1
        for i = 2:N+1
            %disp(i);
            % NOTE, THIS IS BASICALLT GAUSS SEIDEL ITERATION
            %disp(f_c);
            c(i,j) = (D_c*t_size/(h_size)^2) * ( c(i+1,j) + c(i,j+1) - 4*c(i,j) + c(i-1,j) + c(i,j-1) ) + c(i,j) - f_c(i,j)*t_size;
            %disp(c(i,j));
        end
    end
    
    %disp(c_in(1,:));
    % SETS BOUNDARY CONDITIONS
%     c(1,:)=c_in(1,:);
%     c(N,:)=c_in(1,:);
%     c(:,1)=c_in(1,:);
%     c(:,N)=c_in(1,:);
end

