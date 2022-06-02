function grid_tot = grid_total(grid1,grid2)

N = length(grid1(:,1)); % length of grid (number of cells)
grid_tot = -2*ones(N,N);

for i = 2:N-1
    for j = 2:N-1
        grid_tot(i,j) = (grid1(i,j) == 1 | grid1(i,j) == 2 | grid1(i,j) == -1) + (grid2(i,j) == 1 | grid2(i,j) == 2 | grid2(i,j) == -1);
        if or(grid1(i,j) == 3, grid2(i,j) == 3)
            grid_tot(i,j) = 3;
        end
        if grid1(i,j) == 3 && grid2(i,j) == 3
            grid_tot(i,j) = 4;
        end
    end
end


end

