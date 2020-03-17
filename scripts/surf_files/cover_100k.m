% Just plotting

load('re100k_surface.mat');


LoadCmaps
%CmapNames'
cmap = cmaps{1};

%[col]=cbrewer('div', 'RdGy', 5000, 'pchip');
%cmap.colormap=col;

[mapsize, ~]=size(cmap.colormap);

l1 = 1:mapsize;
ncols =5000;
l2 = linspace(1,mapsize,ncols);
r = interp1(l1,cmap.colormap(:,1),l2,'pchip');
g = interp1(l1,cmap.colormap(:,2),l2,'pchip');
b = interp1(l1,cmap.colormap(:,3),l2,'pchip');
cmap.colormap = [r' g' b'];

cmap.colormap = flipud(cmap.colormap);


% Flags/parameters

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
iftrgrey=0;                 % Plot transition location contour
ifsave = 0;
axfs=24;                % axis font size
lafs=36;
figpos = [0.10 0.10 0.70 0.6];

%plotaspect = [1.0 1.7 1.0];

close all

figure(2)
set(gcf, 'Units', 'normalized')
set(gcf, 'OuterPosition', figpos)
set(gcf,'Renderer', 'opengl')

i=0;
if (npts11>0)
  i=i+1;
  splot(i)=surf((surf_t11-ptch_start)/Tosc,surf_x11,surf_v11,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
  view(2)
%  xlabel('$\frac{t}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
%  ylabel('$x/c$', 'FontSize', lafs)
  ylim([0 1])
  hold on

  if (ifcontour)
    cplot{i}=contour((surf_t11-ptch_start)/Tosc,surf_x11,surf_v11, [0 0], 'LineColor', 'k', 'LineWidth', 1.0 );
  end 

end
%colormap(cmap.colormap);

%colorbar;
axis tight
set(gca,'YDir','reverse')
%set(ax1, 'YTickMode', 'manual')
%yticks = [0:20]*0.25;
%set(ax1, 'YTick', yticks);
%set(gca, 'FontSize', axfs)
%set(ax1, 'PlotBoxAspectRatio', [1 1.5 1])
%pbaspect(ax1,plotaspect)

axis off

%f2pos = get(gcf,'OuterPosition');
%axpos=get(ax1,'Position');
%axpos(2)=axpos(2)-0.10;
%axpos(4)=axpos(4)+0.08;
%set(ax1,'Position', axpos);


svfname = ['cover_100k.png'];
if (ifsave)
  SaveFig(gcf, svfname, destn, 1)
end




