function f_c = oxy_rate(R, grid_type)
h_size = 25*10^-4; % room step size in cm (real cell size)
k_r = 6;
T_r = 0.675;
tau = 16*60*60; % cell cycle time in seconds(16 h)
n_0 = h_size^(-2);
c_0 = 1.7*10^-8;

r_c = 2.3*10^-16*tau*n_0/c_0; % base oxygen consumption

%disp("R");
%disp(R);
%disp("TEST");
%disp(grid_type);

%GLU_RATE Summary of this function goes here
%   Detailed explanation goes here
    if grid_type == 2 % If cell is anaerobic
        f_c = 0;
    elseif grid_type == 1 % If cell is aerobic
        F = max(k_r*(R - T_r)+1, 0.25); % Oxygen consumption
        f_c = r_c*F;
    elseif grid_type == 0 % If cell is empty
        f_c = 0;
    elseif grid_type == -1 % If cell is dead
        f_c = 0;
    end
    
    if grid_type == 3 % If cell is healthy
        f_c = max(k_r*(R - T_r)+1, 0.25);
    end
end

