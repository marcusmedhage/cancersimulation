function out = cell_resp(oxy ,glu, pH, nbhr)
c_sap = 0.15; % threshold for apoptosis oxygen
c_gap = 0.5; % threshold for apoptosis glucose
n_max = 8; % threshold for maximum neighbors for proliferation
c_hap = 7.1;

out = zeros(length(oxy), 1);
for n = 1:length(oxy)
    
    if oxy(n) <= c_sap || glu(n) <= c_gap || pH(n) <= c_hap % apoptosis
        val = -1;
    elseif nbhr(n) <= n_max % proliferate, if less neighbors
        val = 1;
    else
        val = 0; % Quiescent
    end
    out(n) = val;
end

