% Just plotting

load('re750k_surface.mat');


% Flags/parameters

%fol = 're750k_pitch';
%ifhdr = 1;
axfs = 9;               % axis fontsize
%lfs = 16;               % legend fontsize
%ifcols = 1;
%ifplot = 0;             % plot individual wall profiles
%tlast = 6.00;          % start from this time
%tstart0 = tlast;
%tend = 100;             % stop at this time
destn = 'plots_test/';
ifcontour=0;            % make contour plot for zero shear stress.
iftr = 1;               % overlay transition points on shear stress space-time plot
iftrportrait=1;         % Plot transition phase portrait
ifxfoil=1;                    % plot xfoil transition location data
ifsave = 1;             % Save space-time plots
ifczplot = 0;           % plot normal force variation
ifczsave = 0;


%ifcp = 0;
lafs=14;                % Latex font size
ifflip=1;
cbloc='Northoutside';
cbheight = 0.75;
%Tosc=1;                 % temporary
%ptch_start=6.0;
%phase=0;


npts8=0;                % remove higher order results

close all

figure(2)
i=0;
%ax1=axes('Position', [0.25 0.1548 0.6589 0.7702]);
ax1=axes;
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
axis tight
hold on
set(ax1, 'YTickMode', 'manual')
yticks = [0:20]*0.5;
set(ax1, 'YTick', yticks);
set(ax1, 'FontSize', axfs)
set(ax1, 'PlotBoxAspectRatio', [1 2.0 1])

if ~ifflip
  ylbl = get(ax1,'YLabel');
  ylblpos = get(ylbl,'Position');
  set(ylbl, 'Position', ylblpos + [-0.01 0 0])
end


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
  set(ylbl,'FontSize', lafs)
  xlabel(ax1,'$x/c$')
  xlbl=get(ax1,'Xlabel');
  set(xlbl,'FontSize', lafs)

  axpos=get(ax1,'Position');
  axpos(2)=axpos(2)-0.10;
  axpos(4)=axpos(4)+0.08;
  set(ax1,'Position', axpos);

  pause(2)   
  SaveFig(gcf, svfname, destn, 1)
end


if (iftr)
  tr = load('tr750k.mat');
  figure(2)
  if npts8>0
    cfmax = max(max(surf_v8));
  else
    cfmax = max(max(surf_v5));
  end
  if ifflip
    cfmax=-cfmax;
  end  

  zdata = zeros(length(tr.tr_time),1)+cfmax;
  trplot = plot3(ax1,tr.trx_uv, (tr.tr_time-ptch_start)/Tosc,zdata, 'LineWidth', 2, 'Color', 'm');
  svfname = ['cf_time_surf750k_tr.eps'];


  if (ifsave)
    if ifflip
%      view([90 -90])
%      set(cb1,'Position',cbpos); 
%      set(cb1,'Location', cbloc); 

      set(ax1,'Position',axpos);
      set(cb1,'Position',cbpos);
    end
 
    SaveFig(gcf, svfname, destn, 1)
  end
end      

%return
%cplot=contour(ax1,surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
%view(2)
%svfname = ['cf_time_surf_contour.eps'];
%SaveFig(gcf, svfname, destn, 1)

figure(3)
ax3=axes;
ax4=axes;
axes(ax3);
j=0;
if (npts5>0)
  axes(ax3)
  j=j+1;
  gplot(j)=surf(ax3,surf_x5,(surf_t5-ptch_start)/Tosc-0.0,surf_c5,surf_c5,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
  set(ax3,'Color', 'none')
  view(2)
  colormap(ax3,'gray');
  %ylabel('t/T_{osc}', 'FontSize', 16)
  xlabel('$x/c$', 'FontSize', lafs)
  xlim([0 1])
  axis tight
  hold on

  axes(ax4)
  gplot2(j)=surf(ax4,surf_x5,(surf_t5-ptch_start)/Tosc,surf_c5,surf_c5,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
  view(2)
  xlim([0 1])
  axis tight
  colormap(ax4,'gray');

end

if (npts8>0)
  axes(ax3)    
  j=j+1;
  gplot(j)=surf(ax3,surf_x8,(surf_t8-ptch_start)/Tosc-0.0,surf_c8,surf_c8,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
  set(ax3,'Color', 'none')
  view(2)
  colormap(ax3,'gray');
  %ylabel('t/T_{osc}', 'FontSize', 16)
  xlabel('$x/c$', 'FontSize', lafs)
  xlim([0 1])
  axis tight
  hold on

  axes(ax4)
  gplot2(j)=surf(ax4,surf_x8,(surf_t8-ptch_start)/Tosc,surf_c8,surf_c8,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
  view(2)
  xlim([0 1])
  axis tight
  colormap(ax4,'gray');
 
end
set(ax3, 'YTickMode', 'manual')
yticks = [0:20]*0.25;
set(ax3, 'YTick', yticks);
set(ax3, 'FontSize', 15)

lncol1 = 'blue';
lncol2 = 'red';

nlines = floor((tlast-tstart0)/Tosc*4);
for i=1:nlines
  if (mod(i,2)==1)
    icol = lncol1;
  else
    icol = lncol2;
  end
  if ifflip
    zpts=-[2 2];
  else
    zpts= [2 2];
  end

  ypts = [i i]*0.25;
  iln(i) = line([0 1], ypts, zpts, 'LineStyle', '--', 'LineWidth', 1.0, 'Color', icol, 'Parent', ax4);

end

set(ax4,'YAxisLocation', 'right');
set(ax4, 'XTick', []);
set(ax4, 'Color', 'none')
set(ax4, 'YTickMode', 'manual')
yticks = [0:nlines]*0.25;
set(ax4, 'YTick', yticks);
ylbl = {'3p/2', '0', 'p/2', 'p'};
set(ax4, 'YTickLabel', ylbl, 'FontName', 'symbol', 'FontSize', axfs);
%ylabel('Oscillation phase', 'Font','Helvetica', 'FontSize', 16)
%set(ax4, 'PlotBoxAspectRatio', [1 1.5 1])

%colorbar
set(ax3,'OuterPosition', get(ax1,'OuterPosition'));
set(ax4,'OuterPosition', get(ax1,'OuterPosition'));

set(ax3,'Position', get(ax1,'Position'));
set(ax4,'Position', get(ax1,'Position'));

figure(3)
svfname = ['cf_time_surf_grey750k.eps'];
if (ifsave)
  if ifflip
    view(ax3,[90 -90])
    view(ax4,[90 -90])
  end    
  SaveFig(gcf, svfname, destn, 1)
end

%% pressure
figure(4)
k=0;
ax1=axes;
if (npts5>0)
  k=k+1;
  pplot(i)=surf(ax1,surf_x5,(surf_t5-ptch_start)/Tosc-0.0,surf_p5,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
  view(2)
  ylabel('$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
  xlabel('$x/c$', 'FontSize', lafs)
  xlim([0 1])
  hold on

end

if (npts8>0)
  k=k+1;
  pplot(k)=surf(ax1,surf_x8,(surf_t8-ptch_start)/Tosc-0.0,surf_p8,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
  view(2)
  ylabel('$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
  xlabel('$x/c$', 'FontSize', lafs)
  xlim([0 1])
  hold on

end
colorbar;
axis tight
hold on
set(ax1, 'YTickMode', 'manual')
yticks = [0:20]*0.25;
set(ax1, 'YTick', yticks);
set(ax1, 'FontSize', axfs)
%set(ax1, 'PlotBoxAspectRatio', [1 1.5 1])
figure(4)
svfname = ['cp_time_surf750k.eps'];
if (ifsave)
  if ifflip
    view([90 -90])
  end    

  SaveFig(gcf, svfname, destn, 1)
end


% cz plot
area=0.15*1;
U0=1.0;
rho=1.0;
norm = 0.5*rho*U0^2*area;
alpha_0 = 3.4;
dalpha  = ptch_amp;    


if ifczplot

  % time series    
  figure(5)
  plot_czt = plot((cz_time-ptch_start)/Tosc,cz/norm, 'LineWidth', 2);
  ylabel('$C_{L}$', 'FontSize', 24)    
  xlabel('$t-t_{0}/T_{osc}$', 'FontSize', 24)    

  svfname = ['cz_time750k.eps'];
  if (ifczsave)
    SaveFig(gcf, svfname, destn, 1)
  end

  % Phase potrait 
  alpha = alpha_0 + dalpha*sin(omega*(cz_time-ptch_start) + phase);
  ind=find(cz_time>31.5);

  figure(6)
  plot_phasepotrait = plot(alpha,cz/norm, 'LineWidth', 2); hold on
  plot_phasepotrait2 = plot(alpha(ind),cz(ind)/norm, 'r', 'LineWidth', 2); hold on
  ylabel('$C_{L}$', 'FontSize', 24)    
  xlabel('$\alpha^{\circ}$', 'FontSize', 24)

  svfname = ['cz_alpha750k.eps'];
  if (ifczsave)
    SaveFig(gcf, svfname, destn, 1)
  end

end

if iftrportrait
  figure(7)
  plot(tr.alpha, tr.trx_uv, 'b', 'LineWidth', 2); hold on
  plot(tr.alpha, tr.trx_ww, 'r', 'LineWidth', 2);
  smoothpvar=tr.trx_ww;
  for jj=1:5
    span=10;
    smoothpvar = smooth(smoothpvar,span);
  end
  ph = arrowh(tr.alpha,smoothpvar,'k',[600,90],[8 15 20 28]);
  xlabel('$\alpha[^{\circ}]$')    
  ylabel('$x/c$')
  set(gca,'FontSize', axfs)
  legend({'$\overline{u''v''}$', '$\overline{w''w''}$'}, 'Interpreter', 'latex','FontSize', 22, 'Location', 'Best')

  set(gca,'FontSize', axfs)
  ylbl=get(gca,'Ylabel');
  set(ylbl,'FontSize', lafs)
  xlbl=get(gca,'Xlabel');
  set(xlbl,'FontSize', lafs)


  svfname = ['750k_transition_alpha.eps'];
  SaveFig(gcf, svfname, destn, 1)

% phase lagged portrait
  phase_lag=-62.89*pi/180;
  tr.alpha2 = alpha_0 + dalpha*sin(omega*(tr.tr_time-ptch_start) + phase + phase_lag);

  figure(8)
  plot(tr.alpha2, tr.trx_uv, 'b', 'LineWidth', 2); hold on
%  plot(tr.alpha2, tr.trx_ww, 'r', 'LineWidth', 2);
  smoothpvar=tr.trx_ww;
  for jj=1:5
    span=10;
    smoothpvar = smooth(smoothpvar,span);
  end
%  ph = arrowh(tr.alpha2,smoothpvar,'k',[600,90],[12 17 29]);
  xlabel('$\alpha_{e}[^{\circ}]$', 'FontSize', lafs)    
  ylabel('$x/c$', 'FontSize', lafs)
  set(gca,'FontSize', axfs)
  ylbl=get(gca,'Ylabel');
  set(ylbl,'FontSize', lafs)
  xlbl=get(gca,'Xlabel');
  set(xlbl,'FontSize', lafs)

%  legend({'$\overline{u''v''}$', '$\overline{w''w''}$'}, 'FontSize', 22, 'Location', 'Best')

  if ifxfoil
    re750 = importdata('test_ed36f128+14_re7.5e5.dat');
    
%    tr_re750_24 = interp1(re750.data(:,1),re750.data(:,7),2.4);
%    tr_re750_44 = interp1(re750.data(:,1),re750.data(:,7),4.4);
    ind1=re750.data(:,1)>=1.5;
    ind2=re750.data(:,1)<5;
    ind3=find(ind1.*ind2);
    
    figure(8)
    plot(re750.data(ind3,1),re750.data(ind3,7), '--k', 'LineWidth', 3); hold on
  end
  svfname = ['750k_transition_alpha_e.eps'];
  SaveFig(gcf, svfname, destn, 1)



end







