function grid = Hard_limits_check(grid,c_disc,g_disc,pH_disc)
% Function that checks for hard limits for oxygen, glucose and pH.

N = length(grid(:,1));

c_ap = 0;
g_0 = 1.3*10^-8; % background concentration glucose
g_ap = 0;%6.5*10^-9/g_0; % apoptosis threshhold for glucose
pH_ap = 6.5;

    for i = 2:N-1
        for j = 2:N-1
            
            if or(grid(i,j) == 1,grid(i,j) == 2) % for aerobic cell
                if pH_disc(i,j) < pH_ap % If pH too high
                     grid(i,j) = 0; % cell dies
                end
                if c_disc(i,j) < c_ap
                    grid(i,j) = -1; % cell dies
                end
                if g_disc(i,j) < g_ap
                    grid(i,j) = -1; % cell dies
                end
            end
            
            if grid(i,j) == 3 % for healthy cell
                if pH_disc(i,j) < pH_ap % If pH too high
                     grid(i,j) = 0; % cell dies
                end
                if c_disc(i,j) < c_ap
                    grid(i,j) = -1; % cell dies
                end
                if g_disc(i,j) < g_ap
                    grid(i,j) = -1; % cell dies
                end
            end
            
        end % END for j
    end % END for i
    
end