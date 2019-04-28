% Just plotting

clear
clc
close all

load('re750k_impulsek04.mat');

surf_v9 = log10(abs(surf_v9)+1.0e-15);
surf_p9 = log10(abs(surf_p9)+1.0e-15);


% Flags/parameters

%fol = 're750k_pitch';
%ifhdr = 1;
tlfs = 16;               % title font size
axfs = 16;               % axis fontsize
%lfs = 16;               % legend fontsize
%ifcols = 1;
%ifplot = 0;             % plot individual wall profiles
%tlast = 6.00;          % start from this time
%tstart0 = tlast;
%tend = 100;             % stop at this time
destn = 'plots_test/';
ifcontour=0;            % make contour plot for zero shear stress.
iftr = 0;               % overlay transition points on shear stress space-time plot
iftrportrait=0;         % Plot transition phase portrait
ifxfoil=0;                    % plot xfoil transition location data
ifsave = 1;             % Save space-time plots
ifczplot = 0;           % plot normal force variation
ifczsave = 0;
ifcf = 0;
ifpressure = 1;

%ifcp = 0;
lafs=16;                % Latex font size
ifflip=0;
cbloc='Eastoutside';
cbheight = 0.75;
%Tosc=1;                 % temporary
%ptch_start=6.0;
%phase=0;


%npts8=0;                % remove higher order results

close all

i=0;
%ax1=axes('Position', [0.25 0.1548 0.6589 0.7702]);

if (ifcf)
  figure(2)
  figpos=[0.25 0.25 0.24 0.50];
  set(gcf, 'Units', 'normalized');
  %set(gcf, 'Renderer', 'Painters')
  set(gcf, 'Position', figpos)
  
  ax1=axes;
  % N=6
  if (npts5>0)
    i=i+1;
    splot(i)=surf(ax1,surf_x5,(surf_t5-ptch_start)/Tosc-0.0,surf_v5,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel(ax1,'$\frac{t-t_{0}}{T_{osc}}$','rot', 0)
    xlabel(ax1,'$x/c$', 'FontSize', lafs)
  
    xlim([0 1])
    hold on
  
    if (ifcontour)
      cplot{i}=contour(ax1,surf_x5,(surf_t5-ptch_start)/Tosc,surf_v5, [0 0], 'LineColor', 'k', 'LineWidth', 1.5  );
    end    
  end
  % N=8
  if (npts8>0)
    i=i+1;
    splot(i)=surf(ax1,surf_x8,(surf_t8-ptch_start)/Tosc-0.0,surf_v8,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel(ax1,'$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
    xlabel(ax1,'$x/c$', 'FontSize', lafs)
    xlim([0 1])
    hold on
  
    if (ifcontour)
      cplot{i}=contour(ax1,surf_x8,(surf_t8-ptch_start)/Tosc,surf_v8, [0 0], 'LineColor', 'k', 'LineWidth', 1.5 );
    end 
  
  end
  % N=9
  if (npts9>0)
    i=i+1;
    splot(i)=surf(ax1,surf_x9,(surf_t9-ptch_start)/Tosc-0.0,surf_v9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel(ax1,'$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
    xlabel(ax1,'$x/c$', 'FontSize', lafs)
    xlim([0 1])
    hold on
  
    if (ifcontour)
      cplot{i}=contour(ax1,surf_x9,(surf_t9-ptch_start)/Tosc,surf_v9, [0 0], 'LineColor', 'k', 'LineWidth', 1.5 );
    end 
  
  end
  
  %set(ax1,'ZScale', 'log')
  
  %axis tight
  hold on
  %set(ax1, 'YTickMode', 'manual')
  %yticks = [0:20]*0.5;
  %set(ax1, 'YTick', yticks);
  %set(ax1, 'FontSize', axfs)
  %set(ax1, 'PlotBoxAspectRatio', [1 2.0 1])
  
  %if ~ifflip
  %  ylbl = get(ax1,'YLabel');
  %  ylblpos = get(ylbl,'Position');
  %  set(ylbl, 'Position', ylblpos + [-0.01 0 0])
  %end
  colormap(jet);
  colorbar
  
  
  svfname = ['cf_time_surf750k.eps'];
  if (ifsave)
    if ifflip
      view([90 -90])
      cb1=colorbar('FontSize', axfs);  
      set(cb1,'Location', cbloc);
      cbpos= get(cb1,'Position');
  %    ax1pos=get(ax1,'Position');   
  %    cbpos(1) = ax1pos(1);
  %%    cbpos(3) = ax1pos(3);
      cbpos(2) = cbpos(2)-0.02;  
      set(cb1,'Position',cbpos);
    else
      cb1=colorbar('FontSize', axfs);  
    end
    ylbl=get(ax1,'Ylabel');
    set(ylbl,'FontSize', lafs+6)
    xlabel(ax1,'$x/c$')
    xlbl=get(ax1,'Xlabel');
    set(xlbl,'FontSize', lafs)
  
    if ifflip    
      axpos=get(ax1,'Position');
      axpos(2)=axpos(2)-0.10;
      axpos(4)=axpos(4)+0.08;
      set(ax1,'Position', axpos);
    end
  
    pause(2)   
    SaveFig(gcf, svfname, destn, 1)
  end
end   % ifcf

%return
%cplot=contour(ax1,surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
%view(2)
%svfname = ['cf_time_surf_contour.eps'];
%SaveFig(gcf, svfname, destn, 1)

%% pressure
if ifpressure
  figure(4)
  figpos=[0.35 0.25 0.24 0.50];
  set(gcf, 'Units', 'normalized');
  %set(gcf, 'Renderer', 'Painters')
  set(gcf, 'Position', figpos)

  k=0;
  ax2=axes;
  if (npts5>0)
    k=k+1;
    pplot(i)=surf(ax2,surf_x5,(surf_t5-ptch_start)/Tosc-0.0,surf_p5,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel('$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs)
    xlabel('$x/c$', 'FontSize', lafs)
    xlim([0 1])
    hold on
  
  end
  
  if (npts8>0)
    k=k+1;
    pplot(k)=surf(ax2,surf_x8,(surf_t8-ptch_start)/Tosc-0.0,surf_p8,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel('$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs)
    xlabel('$x/c$', 'FontSize', lafs)
    xlim([0 1])
    hold on
  
  end
  
  if (npts9>0)
    k=k+1;
    pplot(k)=surf(ax2,surf_x9,(surf_t9-ptch_start)/Tosc-0.0,surf_p9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel('$T$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs)
    xlabel('$x/c$', 'FontSize', lafs)
    xlim([0 1])
    hold on
  
  end
  title('$k=1.0$','FontSize',tlfs)
  
  colormap(jet);
  cb2=colorbar('peer', ax2);
  set(cb2,'Location', cbloc);
  axis tight
  hold on
%  set(ax2, 'YTickMode', 'manual')
%  nphase=2;
%  yticks = [0:30]*1/nphase;
%  set(ax2, 'YTick', yticks);
%  set(ax2, 'FontSize', axfs)
  %set(ax1, 'PlotBoxAspectRatio', [1 1.5 1])
  figure(4)
  svfname = ['cp_time_surfk04.png'];
  pause(2)
  if (ifsave)
    if ifflip
      view([90 -90])
    end    
  
    SaveFig(gcf, svfname, destn, 1)
  end
end         % ifpressure








