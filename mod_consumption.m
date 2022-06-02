function F = mod_consumption(R)
%MOD_CONSUMPTION Summary of this function goes here
%   Detailed explanation goes here
k_r = 6; %modulation strength
T_r = 0.675; % target response;
F = max(k_r*(R - T_r)+1, 0.25);
end

