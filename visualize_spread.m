function [outputArg1,outputArg2] = visualize_spread(c_disc, g_disc, h_disc, pH_disc, grid1, grid2, grid_tot, state1, state2,num_of_plots,t_step)
%VISUALIZE_SPREAD Summary of this function goes here
%   Detailed explanation goes here
    
%     vidfile1 = VideoWriter('Concentrations.avi','Motion JPEG AVI');
%     vidfile1.FrameRate = 5;
%     open(vidfile1);
% 
%     vidfile2 = VideoWriter('Grids.avi','Motion JPEG AVI');
%     vidfile2.FrameRate = 5;
%     open(vidfile2);

    T = size(c_disc,3);

    for k = 1:num_of_plots/t_step:T-1

        cycles = (k-1)*t_step; % number of cycles

%         figure('visible','off')
        figure(1)
        
        subplot(2,2,1)
        surf(c_disc(:,:,k))
        drawnow
        title("Oxygen after " +cycles+ " cycles")
    
        subplot(2,2,2)
        surf(g_disc(:,:,k))
        drawnow
        title("Glucose after " +cycles+ " cycles")
    
        subplot(2,2,3)
        surf(h_disc(:,:,k))
        drawnow
        title("Hydrogen after " +cycles+ " cycles")
    
        subplot(2,2,4)
        surf(pH_disc(:,:,k))
        drawnow
        title("pH after " +cycles+ " cycles")

%         F(k) = getframe(gcf); 
%         writeVideo(vidfile1,F(k));

%         figure('visible','off')
        figure(2)

        subplot(1,2,1)
        imagesc(grid1(:,:,k));
        colormap(subplot(1,2,1),[0 0 0;1 0 0; 1 1 1; 0 1 0; 1 1 0; 0.9608    0.9608    0.8627]);
        set(gca, 'CLim', [-2 3]);
        drawnow
        title("Grid after " +cycles+ " cycles")

        subplot(1,2,2)
        imagesc(state1(:,:,k));
        colormap(subplot(1,2,2),[0 0 0;1 0 1; 0 1 1; 0 0 1; 1 1 1; 0.9608    0.9608    0.8627]);
        set(gca, 'CLim', [-2 3]);
        drawnow
        title("State after " +cycles+ " cycles")
        
        figure(3)

        subplot(1,2,1)
        imagesc(grid2(:,:,k));
        colormap(subplot(1,2,1),[0 0 0;1 0 0; 1 1 1; 0 1 0; 1 1 0; 0.9608    0.9608    0.8627]);
        set(gca, 'CLim', [-2 3]);
        drawnow
        title("Grid after " +cycles+ " cycles")

        subplot(1,2,2)
        imagesc(state2(:,:,k));
        colormap(subplot(1,2,2),[0 0 0;1 0 1; 0 1 1; 0 0 1; 1 1 1; 0.9608 0.9608 0.8627]);
        set(gca, 'CLim', [-2 3]);
        drawnow
        title("State after " +cycles+ " cycles")
        
%         F(k) = getframe(gcf); 
%         writeVideo(vidfile2,F(k));

        figure(4)
        
        subplot(1,1,1)
        imagesc(grid_tot(:,:,k));
        colormap(subplot(1,1,1),[0 0 0;1 1 1; 1 1 1; 0 1 0; 1 1 0; 0.9608 0.9608 0.8627; 0.8500 0.3250 0.0980]);
        set(gca, 'CLim', [-2 4]);
        drawnow
        title("Grid after " +cycles+ " cycles")

        pause(0.1),

    end
    
%     close(vidfile1)
%     close(vidfile2)
       
end



