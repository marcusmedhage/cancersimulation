function h = Diff_hydrogen(h_in, f_h, t_end)
% Function that control the glucose concentration after 1 cell cycle, given
% input hydrogen matrix h_in, number of cells N, input grid and conumption R.

h_size = 25*10^-4; % room step size in cm (real cell size)
tau = 16*60*60; % cell cycle time in seconds(16 h)
t_size = 0.001/tau; % time step size in seconds
T = t_end; % end time (1 cycle)
t = 0:t_size:T; % time discretization
Ts = length(t);

L = 1;
n_0 = h_size^(-2);
h_0 = 1.0*10^-13;
N = size(h_in, 1)-2;


D_h = 1.1*10^-5*tau/L^2; % hydrogen diffusion constant
r_h = 1.5*10^-18*tau*n_0/h_0; % base hydrogen increase per cycle
k_r = 6;
T_r = 0.675;


%Euler forward
h = h_in; % intitial conditions
%disp(Ts);

for k = 2:Ts
    for j = 2:N+1
        for i = 2:N+1
            h(i,j) = (D_h*t_size/(h_size)^2) * ( h(i+1,j) + h(i,j+1) - 4*h(i,j) + h(i-1,j) + h(i,j-1) )+  h(i,j) + f_h(i,j)*t_size;
            %disp(h(i,j));
        end
    end
    
 %       disp(h(1,:));
%     h(1,:)=h(1,:);
%     h(N,:)=h(1,:);
%     h(:,1)=h(1,:);
%     h(:,N)=h(1,:);
end

end
