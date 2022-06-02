function processed = inputProcess(x)
%INPUTPROCESS Summary of this function goes here
%   Detailed explanation goes here

%OXY, GLU, PH, NEIGHBOURS
%
xmin = [0 0 5 0]';
xmax = [1 1 9 9]';
ymax = 1;
ymin = -1;

processed = ((ymax - ymin).*(x - xmin)./(xmax- xmin)) + ymin;

end

