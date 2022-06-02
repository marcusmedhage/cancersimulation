function f_h = hyd_rate(R, grid_type)
h_size = 25*10^-4; % room step size in cm (real cell size)
k_r = 6;
T_r = 0.675;
tau = 16*60*60; % cell cycle time in seconds(16 h)
n_0 = h_size^(-2);
h_0 = 1.0*10^-13;

r_h = 1.5*10^-18*tau*n_0/h_0; % base hydrogen increase per cycle

%GLU_RATE Summary of this function goes here
%   Detailed explanation goes here
    if grid_type == 2 % If cell is alive and anaerobic
        F = max(k_r*(R - T_r)+1, 0.25); % hydrogen increase
        f_h = r_h*F;
    elseif grid_type == 1 % If cell is alive and aerobic
        f_h = 0;
    elseif grid_type== 0 % If cell is empty
        f_h = 0;
    elseif grid_type == -1 % If cell is dead
        f_h = 0;
    end
    
    if grid_type == 3 % If cell is healthy
        f_h = max(k_r*(R - T_r)+1, 0.25);
    end
    
end