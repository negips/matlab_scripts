function  nekchecks(El,Nx,Ny,nelx,nely,nelv,plotgll)

     cmap = lines(nelv);

     h1=figure;               % Plot the grid
     hold on
     for ii=1:nelv
          plot(El(ii).xc([1,2]),El(ii).yc([1,2]), 'Color', cmap(ii,:), 'LineWidth', 2);
          plot(El(ii).xc([2,3]),El(ii).yc([2,3]), 'Color', cmap(ii,:), 'LineWidth', 2);
          plot(El(ii).xc([3,4]),El(ii).yc([3,4]), 'Color', cmap(ii,:), 'LineWidth', 2);
          plot(El(ii).xc([4,1]),El(ii).yc([4,1]), 'Color', cmap(ii,:), 'LineWidth', 2);
     end   
     title('Domain grid');
     xlabel('x');
     ylabel('y');

     filename = 'grid_N4_nelv4_uniform';
     destn = './plots';
%     SaveFig(h1, filename, destn, 0)

%%     Gll points distribution
     if plotgll
          h2 =  figure;       % Gll pts
          hold on;
          for ii=1:nelv
               if mod(ii,2) == 1
                    mrker = '-s';
               else
                    mrker = '-o';
               end
               gll_h = plot(El(ii).xm1,El(ii).ym1, mrker, 'Color', cmap(ii,:), 'MarkerSize', 16, ...
                         'LineWidth', 2.5);
          end
     end
%% Initial solution
     h_uic = figure;
     hold on;
     for ee=1:nelv
          surf(El(ee).xm1,El(ee).ym1,El(ee).un,'EdgeColor', 'none');
          colorbar;
          view([0 90]);
     end



