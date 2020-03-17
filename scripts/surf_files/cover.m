% Just plotting

clear
clc
close all

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/surf_files/cbrewer'

sdata100=load('re100k_surface.mat');
sdata750=load('re750k_surfaceN9.mat');

Re=2;

if Re == 1
  surf_t          = sdata100.surf_t11;
  surf_x          = sdata100.surf_x11;
  surf_v          = sdata100.surf_v11;
  ptch_start      = sdata100.ptch_start;
  Tosc            = sdata100.Tosc;

  svfname = ['cover100k.png'];
else
  surf_t          = sdata750.surf_t9;
  surf_x          = sdata750.surf_x9;
  surf_v          = sdata750.surf_v9;
  ptch_start      = sdata750.ptch_start;
  Tosc            = sdata750.Tosc;
  svfname = ['cover750k.png'];

end  

LoadCmaps
%CmapNames'
cmap = cmaps{1};

[col]=cbrewer('div', 'RdBu', 5000, 'pchip');
cmap.colormap=col;

%cmap.colormap = flipud(jet(5000));

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
i=i+1;
splot(i)=surf((surf_t-ptch_start)/Tosc,surf_x,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
view(2)
 xlabel('$\frac{t}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
 ylabel('$x/c$', 'FontSize', lafs)
ylim([0 1])
hold on

if (ifcontour)
  cplot{i}=contour((surf_t-ptch_start)/Tosc,surf_x,surf_v, [0 0], 'LineColor', 'k', 'LineWidth', 1.0 );
end 

colormap(cmap.colormap);

axis tight
set(gca,'YDir','reverse')

axis off
set(gca,'Position', [0 0 1 1]);


if (ifsave)
  SaveFig(gcf, svfname, destn, 1)
end




