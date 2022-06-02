function g = Diff_glucose(g_in, f_g, t_end)
% Function that control the glucose concentration after 1 cell cycle, given
% input glucose matrix g_in, number of cells N, input grid and conumption R.

h_size = 25*10^-4; % room step size in cm (real cell size)
tau = 16*60*60; % cell cycle time in seconds(16 h)
t_size = 0.001/tau; % time step size in seconds
T = t_end; % end time (1 cycle)
t = 0:t_size:T; % time discretization
Ts = length(t);

L = 1;
n_0 = h_size^(-2);
g_0 = 1.3*10^-8;
N = size(g_in, 1)-2;


D_g = 9.1*10^-5*tau/L^2; % glucose diffusion constant
r_g_a = 3.8*10^-17*tau*n_0/g_0; % base glucose consuption (aerobic)
r_g_an = 6.9*10^-16*tau*n_0/g_0; % base glucose consuption (anaerobic)
k_r = 6;
T_r = 0.675;

%Euler forward
g = g_in; % intitial conditions

for k = 2:Ts
    for j = 2:N+1
        for i = 2:N+1
            %disp(g(i,j));
            g(i,j) = (D_g*t_size/(h_size)^2) * ( g(i+1,j) + g(i,j+1) - 4*g(i,j) + g(i-1,j) + g(i,j-1) ) +g(i,j) - f_g(i,j)*t_size;
            %disp(g(i,j));
            %disp(g(i,j));
        end
    end
%     g(1,:)=g(1,:);
%     g(N,:)=g(1,:);
%     g(:,1)=g(1,:);
%     g(:,N)=g(1,:);
end
