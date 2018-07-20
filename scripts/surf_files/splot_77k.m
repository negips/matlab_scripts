% Just plotting


kred=0.0477;
load(['re77k_k' num2str(kred*10000) '_surface.mat']);

% Flags/parameters
U0=1.;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=2*pi/omega;
ptch_amp = 5.5;
ptch_start = 16.0;
axis_x0 = 0.186;
axis_y0 = 0.0;
phase=0;



% fs = 16;                % fontsize
% lfs = 16;               % legend fontsize
% ifcols = 1;
% ifplot = 0;
% tlast = 19.05;
% tstart0 = tlast; 
destn = 'plots_test/';
ifcontour=1;
iftr=0;                 % Plot transition location contour
% iftrabs=0;              % plot point of absolute instability
iftrgrey=1;                 % Plot transition location contour
ifsave = 1;
axfs=24;                % axis font size
lafs=36;
figpos = [0.10 0.10 0.30 0.8];

plotaspect = [1.0 1.7 1.0];

close all

figure(2)
set(gcf, 'Units', 'normalized')
set(gcf, 'OuterPosition', figpos)
set(gcf,'Renderer', 'opengl')

i=0;
ax1=axes;
if (npts7>0)
  i=i+1;
  splot(i)=surf(ax1,surf_x7,(surf_t7-ptch_start)/Tosc-0.0,surf_v7,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
  view(2)
  ylabel(ax1,'$\frac{t}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
  xlabel(ax1,'$x/c$', 'FontSize', lafs)
  xlim([0 1])
  hold on

  if (ifcontour)
    cplot{i}=contour(ax1,surf_x7,(surf_t7-ptch_start)/Tosc,surf_v7, [0 0], 'LineColor', 'k', 'LineWidth', 1.5  );
  end    
end

if (npts9>0)
  i=i+1;
  splot(i)=surf(ax1,surf_x9,(surf_t9-ptch_start)/Tosc-0.0,surf_v9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
  view(2)
  ylabel(ax1,'$\frac{t}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
  xlabel(ax1,'$x/c$', 'FontSize', lafs)
  xlim([0 1])
  hold on

  if (ifcontour)
    cplot{i}=contour(ax1,surf_x9,(surf_t9-ptch_start)/Tosc,surf_v9, [0 0], 'LineColor', 'k', 'LineWidth', 1.5 );
  end 

end
colorbar;
axis tight
hold on
set(ax1, 'YTickMode', 'manual')
yticks = [0:100]*1.00;
set(ax1, 'YTick', yticks);
set(ax1, 'FontSize', axfs)
%set(ax1, 'PlotBoxAspectRatio', [1 1.5 1])
pbaspect(ax1,plotaspect)

ylbl = get(ax1,'YLabel');
ylblpos = get(ylbl,'Position');
set(ylbl, 'Position', ylblpos + [-0.02 0 0])

ylbl=get(ax1,'Ylabel');
set(ylbl,'FontSize', lafs+12)
xlabel(ax1,'$x/c$')
xlbl=get(ax1,'Xlabel');
set(xlbl,'FontSize', lafs)

f2pos = get(gcf,'OuterPosition');
%axpos=get(ax1,'Position');
%axpos(2)=axpos(2)-0.10;
%axpos(4)=axpos(4)+0.08;
%set(ax1,'Position', axpos);


svfname = ['cf_time_surf77k_k' num2str(kred) '.eps'];
%destn = 'plots/';   
if (ifsave)
  SaveFig(gcf, svfname, destn, 1)
end

%if (iftr)
%  tr = load('tr100k.mat');
%  figure(2)
%  cfmax = max(max(surf_v9));    
%  zdata = zeros(length(tr.tr_time),1)+cfmax;
%  trplot = plot3(ax1, tr.trx_uv, (tr.tr_time-ptch_start)/Tosc,zdata, 'LineWidth', 1.5, 'Color', 'm');
%  svfname = ['cf_time_surf100k_tr.eps'];
%
%  if (iftrabs)
%    tabs = 20.65; 
%    trabs = interp1(tr.tr_time,tr.trx_uv,20.65);
%    absplot = plot3(ax1,trabs,tabs/Tosc, cfmax, 'o', 'MarkerSize', 20, 'LineWidth', 2.5, 'MarkerFaceColor', 'b');
%    svfname = ['cf_time_surf100k_trabs.eps'];
%  end 
%
%  pbaspect(ax1,plotaspect)
%  
%%  destn = 'plots/';   
%  if (ifsave)
%    SaveFig(gcf, svfname, destn, 1)
%  end
%end      



%cplot=contour(ax1,surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
%view(2)
%svfname = ['cf_time_surf_contour.eps'];
%destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)

figure(3)
set(gcf, 'Units', 'normalized')
set(gcf, 'OuterPosition', figpos)
set(gcf,'Renderer', 'opengl')
ax3=axes;
%ax4=axes;
axes(ax3);
j=0;
if (npts7>0)
  axes(ax3)
  j=j+1;
  gplot(j)=surf(ax3,surf_x7,(surf_t7-ptch_start)/Tosc-0.0,surf_c7,surf_c7,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
  set(ax3,'Color', 'none')
  view(2)
  colormap(ax3,'gray');
  %ylabel('t/T_{osc}', 'FontSize', 16)
  xlabel('$x/c$', 'FontSize', lafs)
  xlim([0 1])
  axis tight
  hold on

%  axes(ax4)
%  gplot2(j)=surf(ax4,surf_x7,(surf_t7-ptch_start)/Tosc,surf_c7,surf_c7,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
%  view(2)
%  axis tight
end

if (npts9>0)
  axes(ax3)    
  j=j+1;
  gplot(j)=surf(ax3,surf_x9,(surf_t9-ptch_start)/Tosc-0.0,surf_c9,surf_c9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
  set(ax3,'Color', 'none')
  view(2)
  colormap(ax3,'gray');
  %ylabel('t/T_{osc}', 'FontSize', 16)
  xlabel('$x/c$', 'FontSize', lafs)
  xlim([0 1])
  axis tight
  hold on

%  axes(ax4)
%  gplot2(j)=surf(ax4,surf_x9,(surf_t9-ptch_start)/Tosc,surf_c9,surf_c9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
%  view(2)
%  axis tight
end
set(ax3, 'YTickMode', 'manual')
yticks = [0:100]*1.00;
set(ax3, 'YTick', yticks);
set(ax3, 'FontSize', axfs)

lncol1 = 'blue';
lncol2 = 'red';

%nlines = floor((tlast-tstart0)/Tosc)+1;
%for i=1:nlines
%  if (mod(i,2)==1)
%    icol = lncol1;
%  else
%    icol = lncol2;
%  end
%
%  ypts = 3 + [i i]*0.25;
%  iln(i) = line([0 1], ypts, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', icol, 'Parent', ax4);
%
%end

%set(ax4,'YAxisLocation', 'right');
%set(ax4, 'XTick', []);
%set(ax4, 'Color', 'none')
%set(ax4, 'YTickMode', 'manual')
%yticks = 3 + [0:nlines]*0.25;
%set(ax4, 'YTick', yticks);
%%ylbl = {'0', 'p/2', 'p', '3p/2'};
%%set(ax4, 'YTickLabel', ylbl, 'FontName', 'symbol', 'FontSize', axfs);
%%ylbl = {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$'};
%%set(ax4, 'YTickLabel', ylbl, 'FontName', 'Latex', 'FontSize', axfs);
%
%%ylabel('Oscillation phase', 'Font','Helvetica', 'FontSize', 16)
%%set(ax4, 'PlotBoxAspectRatio', [1 1.5 1])
%
%ylbl = {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$'};
%set(ax4, 'YTickLabel', ylbl, 'TickLabelInterpreter', 'latex', 'FontSize', axfs+4);


%colorbar
%colormap(ax4,'gray');

%set(ax3,'OuterPosition', get(ax1,'OuterPosition'));
%set(ax4,'OuterPosition', get(ax1,'OuterPosition'));

ax1pos=get(ax1,'Position');

%set(ax3,'Position', ax1pos+[0.15 0 -0.29 0]);
%set(ax4,'Position', ax1pos+[0.15 0 -0.29 0]);
set(gcf,'OuterPosition',figpos)
pbaspect(ax3,plotaspect)
%pbaspect(ax4,plotaspect)

svfname = ['cf_time_surf_grey77k_k' num2str(kred) '.eps'];
figure(3)
%destn = 'plots/';   
if (ifsave)
  SaveFig(gcf, svfname, destn, 1)
end

%if (iftrgrey)
%  tr = load('tr100k.mat');
%  figure(2)
%  zdata = zeros(length(tr.tr_time),1)+1.1;
%  trplot2 = plot3(ax4, tr.trx_uv, (tr.tr_time-ptch_start)/Tosc,zdata, 'LineWidth', 1.5, 'Color', 'm');
%  svfname = ['cf_time_surf_grey100k_tr.eps'];
%
%% Remove horizontal lines          
%  for i=1:nlines    
%    set(iln(i),'Visible','off')
%  end    
%  figure(3)
%%  destn = 'plots/';   
%  if (ifsave)
%    SaveFig(gcf, svfname, destn, 1)
%  end
%end      

