function [p,p_e] = pressure_eq(grid1, grid2)
%PRESSURE_EQ Summary of this function goes here
%   Detailed explanation goes here

N = length(grid1);

s_e = (grid1 == 3) + (grid2 == 3); % external pressure sources, where there are healthy cells
s = (grid1 == 1 | grid1 == 2) + (grid2 ==1  | grid2 == 2); % internal pressure, where there are cancer cells

h = 25*10^-4; % room step size in cm (real cell size)

A = zeros(N,N);

for i = 1:N-1
   A(i,i+1) = -1;
   A(i+1,i) = -1;
   A(i,i) = 2;
end
A(N,N) = 2;
%disp(s);
%disp(A);

A = A/h^2;

p = -A\s;

p_e = -A\s_e;

% s = 1 if u = 1, note u=2 in report, but we try u=1 here.
% s = 0 if u = otherwise


end
