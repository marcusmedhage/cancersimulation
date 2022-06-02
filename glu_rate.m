function f_g = glu_rate(R, grid_type)

h_size = 25*10^-4; % room step size in cm (real cell size)
k_r = 6;
T_r = 0.675;
tau = 16*60*60; % cell cycle time in seconds(16 h)
n_0 = h_size^(-2);
g_0 = 1.3*10^-8;

r_g_a = 3.8*10^-17*tau*n_0/g_0; % base glucose consuption (aerobic)
r_g_an = 6.9*10^-16*tau*n_0/g_0; % base glucose consuption (anaerobic)

%OXY_RATE Summary of this function goes here
%   Detailed explanation goes here
    if grid_type == 2 % If cell is alive and anaerobic
        F = max(k_r*    (R-T_r)+1, 0.25); % Glucose consumption
        f_g = r_g_an*F;
        %disp("AN");
    elseif grid_type == 1 % If cell is alive and aerobic
        F = max(k_r*    (R-T_r)+1, 0.25); % Glucose consumption
        f_g = r_g_a*F;
        %disp("AER");
        %disp(F);
        %disp(f_g);
    elseif grid_type == 0 % If cell is empty
        f_g = 0;
    elseif grid_type == -1 % If cell is dead
        f_g = 0;
    end
    
    if grid_type == 3 % If cell is healthy
        f_g = max(k_r*(R - T_r)+1, 0.25);
    end
end

