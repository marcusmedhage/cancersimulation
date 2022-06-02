function grid_out = Dead_cell_removal(grid_in,t_n,k)
% Function that removes dead cells from grid affter t_n cycles.

N = length(grid_in(:,1,1));
k = length(grid_in(1,1,:));

grid_out = grid_in(:,:,end);

 for i = 1:N
     for j = 1:N
            
         if k > t_n % makes sure dead cells are removed after t_n cycles
             dead_cells = 0; % counts number of timesteps been dead
                
             %% POssible to replace with sum(grid[... ...])
             for h = 0:t_n
                 if grid_in(i,j,k-h) == -1
                 dead_cells = dead_cells + 1; % if it was dead earlier
                 end               
             end
             if dead_cells == t_n +1
                 grid_out(i,j) = 0;
             end  % END IF DEAD CELL
         end % END IF k > t_n
            
     end 
 end    