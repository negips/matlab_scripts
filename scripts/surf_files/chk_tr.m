% Checking Transition

clear
clc
close all

load('re750k_surfaceN9.mat');


axfs = 16;               % axis fontsize
%lfs = 16;               % legend fontsize
%ifcols = 1;
%ifplot = 0;             % plot individual wall profiles
%tlast = 6.00;          % start from this time
%tstart0 = tlast;
%tend = 100;             % stop at this time
destn = 'plots_test/';
ifcontour=1;            % make contour plot for zero shear stress.
iftr = 1;               % overlay transition points on shear stress space-time plot
iftrportrait=0;         % Plot transition phase portrait
ifxfoil=0;                    % plot xfoil transition location data
ifsave = 0;             % Save space-time plots
ifczplot = 0;           % plot normal force variation
ifczsave = 0;


%ifcp = 0;
lafs=30;                % Latex font size
ifflip=0;
cbloc='Northoutside';
cbheight = 0.75;
%Tosc=1;                 % temporary
%ptch_start=6.0;
%phase=0;


x2=surf_x9;
y2=surf_t9;
z2=surf_v9;

[r,c]=size(x2);

ind = find(x2<0.25);
z2(ind)=0.;

tol=0.33;
for i=1:r
  [cmax,ind]=sort(z2(i,:),'descend');
  cmax2 = z2(i,ind(1));
  ind2  = ind(1);
  ztmp  = z2(i,1:ind2);
  xtmp  = x2(i,1:ind2);
  ytmp  = y2(i,1:ind2);
  ind3  = find(ztmp<tol*cmax2,1,'last');
  xtr(i) = xtmp(ind3);
  ttr(i) = ytmp(ind3);
  ztr(i) = ztmp(ind3);
  
%  ind3  = find(z2(i,:)>tol*cmax,1,'first');
%  xtr(i) = x2(i,ind3);
%  ttr(i) = y2(i,ind3);
%  ztr(i) = z2(i,ind3);
end

%plot(ttr,xtr)

i=0;
outpos1 = [0.2 0.1 0.3 0.7];
% N=9
if (npts9>0)
  figure(1)
  set(gcf,'Units','Normalized');
  set(gcf,'OuterPosition',outpos1)
  i=i+1;

  splot(i)=surf(surf_x9,(surf_t9-ptch_start)/Tosc-0.0,surf_v9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
  ax1=gca;
  ylabel(ax1,'$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 90, 'FontSize', lafs+6)
  xlabel(ax1,'$x/c$', 'FontSize', lafs)
  view(2)
  xlim([0 1])
  hold on

%  if (ifcontour)
%    cplot{i}=contour(ax1,surf_x9,(surf_t9-ptch_start)/Tosc,surf_v9, [0 0], 'LineColor', 'k', 'LineWidth', 1.5 );
%  end 

end
axis tight
hold on
set(ax1, 'YTickMode', 'manual')
yticks = [0:20]*0.5;
set(ax1, 'YTick', yticks);
set(ax1, 'FontSize', axfs)

if (iftr)
  disp('Transition')
  plot3(xtr,(ttr-ptch_start)/Tosc,ztr+0.001, 'm', 'LineWidth',2);
end  







