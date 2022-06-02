function resp = healthy_resp(oxy, glu, pH)
%HEALTHY_RESP Summary of this function goes here
%   Detailed explanation goes here
    oxy_l = 0.15;
    glu_l = 0.5;
    pH_l = 7.1;

    if oxy > oxy_l && glu > glu_l && pH > pH_l % CELL LIVES AS NORMAL
        resp = 0; % stays quiescent
    else % CELL DIES
        resp = -1;
    end

end

