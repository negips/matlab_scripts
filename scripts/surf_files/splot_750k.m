% Just plotting

clear
clc
close all

load('re750k_surfaceN9.mat');

addpath '~/workstation/git_kth/matlabscripts/scripts/surf_files/cbrewer'

LoadCmaps
%CmapNames'
cmap = cmaps{1};


cover_plot = 1;
if (cover_plot)

  s = surf_v9(:);
  sm = mean(s);
  st = std(s);
  ind = find(surf_v9<(sm-2*st));
  surf_v9(ind)=sm-2*st;
      
  [col]=cbrewer('div', 'RdGy', 5000, 'pchip');
  %cmap.colormap=col;
  
  npt1 = 1500;
  [col1]=cbrewer('div', 'RdBu', npt1, 'pchip');
  
  npt2 = 5000;
  %[col2]=cbrewer('div', 'PuOr', npt2, 'pchip');
  [col2]=cbrewer('div', 'BrBG', npt2, 'pchip');
  
  col = [col2(1:npt2/2,:); col1((npt1/2+1):npt1,:)];
  
  cmap.colormap=col;

end  


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

fontfac=0.5;

%fol = 're750k_pitch';
%ifhdr = 1;
axfs = 28*fontfac;               % axis fontsize
%lfs = 16;               % legend fontsize
%ifcols = 1;
%ifplot = 0;             % plot individual wall profiles
%tlast = 6.00;          % start from this time
%tstart0 = tlast;
%tend = 100;             % stop at this time
destn = 'plots_test/';
%destn = '~/workstation/git_kth/forced_pitching/paper/imgs2/';
ifcontour=1;            % make contour plot for zero shear stress.
iftr = 1;               % overlay transition points on shear stress space-time plot
  trfile='tr750k_n9_2.mat';
<<<<<<< HEAD
=======
iftrportrait=0;         % Plot transition phase portrait
ifxfoil=0;                    % plot xfoil transition location data

ifseparated=0;
ifpressure=0;

>>>>>>> aeeda89880074971def49116caf3607c6d079210
ifsave = 0;             % Save space-time plots
  
iftrportrait=1;         % Plot transition phase portrait
  ifxfoil=1;                          % plot xfoil transition location data
    xfile='saab750k_N8.5.t1.dat';     % file for Xfoil Data
      ifportsave = 1;             % Save space-time plots

ifseparated=0;
ifpressure=0;

ifczplot = 0;           % plot normal force variation
ifczsave = 0;



%ifcp = 0;
<<<<<<< HEAD
lafs=42*fontfac;                % Latex font size
=======
lafs=30;                % Latex font size
>>>>>>> aeeda89880074971def49116caf3607c6d079210
ifflip=1;
cbloc='Northoutside';
cbheight = 0.75;
%Tosc=1;                 % temporary
%ptch_start=6.0;
%phase=0;

%npts8=0;                % remove higher order results

close all

if (~ifflip)
  figpos=[0.25 0.25 0.24 0.50];
else
  figpos=[0.15 0.4 0.65 0.60];
end  

figure(2)
set(gcf, 'Units', 'normalized');
%set(gcf, 'Renderer', 'Painters')
set(gcf, 'OuterPosition', figpos)

i=0;
%ax1=axes('Position', [0.25 0.1548 0.6589 0.7702]);
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
  if (ifflip)
    splot(i)=surf(ax1,(surf_t9-ptch_start)/Tosc-0.0,surf_x9,surf_v9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    xlabel(ax1,'$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
    ylabel(ax1,'$x/c$', 'FontSize', lafs)
    ylim([0 1])
    hold on

    if (ifcontour)
      cplot{i}=contour(ax1,(surf_t9-ptch_start)/Tosc,surf_x9,surf_v9, [0 0], 'LineColor', 'k', 'LineWidth', 1.5 );
    end 

  else            
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

end

axis tight
hold on
if ~(ifflip)
  set(ax1, 'YTickMode', 'manual')
  yticks = [0:20]*0.5;
  set(ax1, 'YTick', yticks);
  set(ax1, 'FontSize', axfs)
else
  set(ax1,'YDir','reverse')
  set(ax1, 'XTickMode', 'manual')
  xticks = [0:40]*0.25;
  set(ax1, 'XTick', xticks);
  set(ax1, 'FontSize', axfs);
  %set(ax1, 'PlotBoxAspectRatio', [1 2.0 1])

  set(ax1, 'YTickMode', 'manual')
  yticks = [0:5]*0.2;
  set(ax1, 'YTick', yticks);
  set(ax1, 'FontSize', axfs)
end  

if ~ifflip
  ylbl = get(ax1,'YLabel');
  ylblpos = get(ylbl,'Position');
  set(ylbl, 'Position', ylblpos + [-0.01 0 0])
end
%colormap(jet);
colormap(cmap.colormap);


svfname = ['cf_time_surf750kN9.png'];
if (ifsave)
  if ifflip
    cb1=colorbar('FontSize', axfs);  
    set(cb1,'Location', cbloc);
    cbpos= get(cb1,'Position');
    ax1pos=get(ax1,'Position'); 
    cbpos(1) = ax1pos(1);
%    cbpos(3) = ax1pos(3);
    cbpos(2) = cbpos(2)+0.01;  
    set(cb1,'Position',cbpos);
    set(cb1,'FontSize',12);

    ylbl=get(ax1,'Ylabel');
    set(ylbl,'FontSize', lafs)
    xlbl=get(ax1,'Xlabel');
    set(xlbl,'FontSize', lafs+6)
   
  else
    cb1=colorbar('FontSize', axfs);
    ylbl=get(ax1,'Ylabel');
    set(ylbl,'FontSize', lafs+6)
    xlbl=get(ax1,'Xlabel');
    set(xlbl,'FontSize', lafs)
  end

  if ifflip    
    axpos=get(ax1,'Position');
    axpos(2)=axpos(2)+0.02;
    axpos(4)=axpos(4)-0.1;
    set(ax1,'Position', axpos);
  end

  pause(2)
  if ~(iftr)
    SaveFig(gcf, svfname, destn, 1)
  end  
end


if (iftr)
  tr = load(trfile);
  figure(2)
  if npts9>0
    cfmax = max(max(surf_v9));
  elseif npts8>0
    cfmax = max(max(surf_v8));
  else
    cfmax = max(max(surf_v5));
  end
%  if ifflip
%    cfmax=-cfmax;
%  end  

  zdata = zeros(length(tr.tr_time),1)+cfmax;
  if ~(ifflip)
    trplot = plot3(ax1,tr.trx_uv,(tr.tr_time-ptch_start)/Tosc,zdata, 'LineWidth', 2, 'Color', 'r');
  else
    trplot = plot3(ax1,(tr.tr_time-ptch_start)/Tosc,tr.trx_uv,zdata, 'LineWidth', 2, 'Color', 'r');
  end      
  svfname = ['cf_time_surf750k_trN9.png'];

  if (ifsave)
    SaveFig(gcf, svfname, destn, 1)
  end
end      

%return
%cplot=contour(ax1,surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
%view(2)
%svfname = ['cf_time_surf_contour.eps'];
%SaveFig(gcf, svfname, destn, 1)


if (ifseparated)
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
  
    axes(ax4)
    gplot2(j)=surf(ax4,surf_x9,(surf_t9-ptch_start)/Tosc,surf_c9,surf_c9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
    view(2)
    xlim([0 1])
    axis tight
    colormap(ax4,'gray');
   
  end
  
  
  set(ax3, 'YTickMode', 'manual')
  nphase=2;
  yticks = [0:100]*1/nphase;
  set(ax3, 'YTick', yticks);
  set(ax3, 'FontSize', 15)
  
  lncol1 = 'blue';
  lncol2 = 'red';
  
  nlines = floor((tlast-ptch_start)/Tosc*nphase);
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
  
    ypts = [i i]*1/nphase;
    iln(i) = line([0 1], ypts, zpts, 'LineStyle', '--', 'LineWidth', 1.0, 'Color', icol, 'Parent', ax4);
  
  end
  
  set(ax4,'YAxisLocation', 'right');
  set(ax4, 'XTick', []);
  set(ax4, 'Color', 'none')
  set(ax4, 'YTickMode', 'manual')
  yticks = [0:nlines]*1/nphase;
  set(ax4, 'YTick', yticks);
  %ylbl = {'3p/2', '0', 'p/2', 'p'};
  %set(ax4, 'YTickLabel', ylbl, 'FontName', 'symbol', 'FontSize', axfs);
  %ylabel('Oscillation phase', 'Font','Helvetica', 'FontSize', 16)
  %set(ax4, 'PlotBoxAspectRatio', [1 1.5 1])
  
  ylbl = {'$3\pi/2$', '$0$', '$\pi/2$', '$\pi$'};
  set(ax4, 'YTickLabel', ylbl,'FontSize', axfs);
  
  
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
end         % ifsepratated




%% pressure
if (ifpressure)
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
  
  if (npts9>0)
    k=k+1;
    pplot(k)=surf(ax1,surf_x9,(surf_t9-ptch_start)/Tosc-0.0,surf_p9,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
    view(2)
    ylabel('$\frac{t-t_{0}}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', lafs+6)
    xlabel('$x/c$', 'FontSize', lafs)
    xlim([0 1])
    hold on
  
  end
  
  colormap(jet);
  colorbar;
  axis tight
  hold on
  set(ax1, 'YTickMode', 'manual')
  nphase=2;
  yticks = [0:30]*1/nphase;
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
end         % ifpressure

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
 
  tr = load(trfile);
     
  gr=(1+sqrt(5))/2;
  figpos=[0.25 0.25 0.35 gr*0.35];
  figure(7)
  set(gcf, 'Units', 'normalized');
  set(gcf, 'OuterPosition', figpos)
  set(gcf, 'Renderer', 'Painters')
%  ax1=axes;
%  axpos=[0.1400 0.1500 0.7650 0.7950];
%  set(ax1,'Position',axpos)
%  set(ax1,'FontSize',axfs)

  plot(tr.alpha, tr.trx_uv, 'LineWidth', 1.5); hold on
%  plot(tr.alpha, tr.trx_ww, 'r', 'LineWidth', 1.5);
  smoothpvar=tr.trx_uv;
  for jj=1:5
    span=10;
    smoothpvar = smooth(smoothpvar,span);
  end
  ph = arrowh(tr.alpha,smoothpvar,'r',[400,90],[8 18 35]);
  xlabel('$\alpha[^{\circ}]$')    
  ylabel('$x_{tr}/c$')
  set(gca,'FontSize', axfs)
  xlim([2.2 4.6])
%  legend({'$\overline{u''v''}$', '$\overline{w''w''}$'}, 'Interpreter', 'latex','FontSize', 22, 'Location', 'Best')

  set(gca,'FontSize', axfs)
  ylbl=get(gca,'Ylabel');
  set(ylbl,'FontSize', lafs)
  xlbl=get(gca,'Xlabel');
  set(xlbl,'FontSize', lafs)

  pause(1)

  if (ifportsave)
    svfname = ['750k_transition_alpha.eps'];
    SaveFig(gcf, svfname, destn, 1)
  end        

% phase lagged portrait
%  phase_lag=-62.89*pi/180;
  phase_lag=-59.00*pi/180;

  tfine = linspace(tr.tr_time(1),tr.tr_time(end),5000);
  vfine = interp1(tr.tr_time,tr.trx_uv,tfine,'pchip');
  ar = areafcn(tfine,vfine,0);
  
  lag0 = -60*pi/180;
  [oplag area exitflag output] = fminsearch(@(lag) areafcn(tfine,vfine,lag),lag0);
%  phase_lag=oplag;
  phase_lag=-61*pi/180;

  tr.alpha2 = alpha_0 + dalpha*sin(omega*(tr.tr_time-ptch_start) + phase + phase_lag);

% Empirically obtained from stats  
  tr24 = 0.766;
  tr34 = 0.659;
  tr44 = 0.436;

  figure(8)
  set(gcf, 'Units', 'normalized');
  set(gcf, 'OuterPosition', figpos)
  set(gcf, 'Renderer', 'Painters')
  plot(tr.alpha2, tr.trx_uv, 'LineWidth', 1.5); hold on
  plot([2.4 3.4 4.4],[tr24 tr34 tr44], 'or', 'MarkerSize', 12, 'LineWidth',3)
  xlabel('$\alpha_{e}[^{\circ}]$', 'FontSize', lafs)    
  ylabel('$x_{tr}/c$', 'FontSize', lafs)
  set(gca,'FontSize', axfs)
  ylbl=get(gca,'Ylabel');
  set(ylbl,'FontSize', lafs)
  xlbl=get(gca,'Xlabel');
  set(xlbl,'FontSize', lafs)

%  legend({'$\overline{u''v''}$', '$\overline{w''w''}$'}, 'FontSize', 22, 'Location', 'Best')

  if ifxfoil
%    re750 = importdata('test_ed36f128+14_re7.5e5.dat');
    re750 = importdata(xfile);
   
%    tr_re750_24 = interp1(re750.data(:,1),re750.data(:,7),2.4);
%    tr_re750_44 = interp1(re750.data(:,1),re750.data(:,7),4.4);
    ind1=re750.data(:,1)>=1.5;
    ind2=re750.data(:,1)<5;
    ind3=find(ind1.*ind2);
    
    figure(8)
    plot(re750.data(ind3,1),re750.data(ind3,7), '--k', 'LineWidth', 3); hold on
  end

  if (ifportsave)
    svfname = ['750k_transition_alpha_e.eps'];
    SaveFig(gcf, svfname, destn, 1)
  end  

end







