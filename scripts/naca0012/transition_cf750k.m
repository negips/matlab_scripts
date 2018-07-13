% Plot cp files

clear
clc
close all

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 're750k_pitch';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
ifplot = 0;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 0.0;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

% Load mean aoa normals
snormals = importdata('surf_normals.24.N5');

x_imp = snormals.data(:,1);
[x_imp I] = sort(x_imp);
y_imp = snormals.data(I,2);

snx = -snormals.data(I,4);
sny = -snormals.data(I,5);

ind = sny>0;
snx_top = snx(find(ind));
sny_top = sny(find(ind));
xt_imp  = x_imp(find(ind));
yt_imp  = y_imp(find(ind));

snx_bot = snx(find(~ind));
sny_bot = sny(find(~ind));
xb_imp  = x_imp(find(~ind));
yb_imp  = y_imp(find(~ind));

for i=1:length(snx_top)
  stx_top(i) = sny_top(i);
  sty_top(i) = -snx_top(i);
end

for i=1:length(snx_bot)
  stx_bot(i) = -sny_bot(i);
  sty_bot(i) = snx_bot(i);
end
%% 

U0=1.;
kred=0.4;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=2*pi/omega;
ptch_amp = 1.0;
ptch_start = 6.;
axis_x0 = 0.35;
axis_y0 = 0.034;
phase=-pi/2;

nplots = 0;
bcnt = 0;
surf_x = [];
surf_t = [];
surf_v = [];
surf_c = [];
xtr = [];
icalld = 0;
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    xmax = max(x(:));
    xmin = min(x(:));      
    Chord=xmax-xmin;

%   trying to find top and bottom elements        
    [val ind] = max(sdata(1).data(:,:,1));
    n12=find(ind==lx1(1));
    n01=find(ind==1);

    y12 = sdata(2).data(:,n12,1);
    y12mean = mean(y12(:));  
    y01 = sdata(2).data(:,n01,1);
    y01mean = mean(y01(:));

    if y12mean>y01mean
      t_els=n12;
      b_els=n01;
    else  
      t_els=n01;
      b_els=n12;
    end

    xtels_mean = mean(sdata(1).data(:,t_els,1));
    [val ind] = max(xtels_mean);
    t_els(ind)=[];                  % remove trailing edge blunt side 

%   Remove a few leading edge elements  
    if (icalld==0)

      for i=1:5   % 5 leading elements
        xtels_mean = mean(sdata(1).data(:,t_els,1));
        [val ind] = min(xtels_mean);
        t_els(ind)=[];                  % remove leading edge element
      end 
    
      total_els = length(t_els);
      icalld = icalld+1;
    else
      ntop = length(t_els);
      while ntop>total_els
        xtels_mean = mean(sdata(1).data(:,t_els,1));
        [val ind] = min(xtels_mean);
        t_els(ind)=[];                  % remove leading edge element 
        ntop = length(t_els);
      end
    end    

      
    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast)
         if ifplot && nplots>0
           delete(pvar)
           delete(trplot) 
         end   
         dtmpx = sdata(1).data(:,t_els,it);
         dtmpy = sdata(2).data(:,t_els,it);
%         dtmp_v = sdata(3).data(:,t_els,it);    
         dtmp_v = -sdata(5).data(:,t_els,it);
         dtmp_cfx = -sdata(5).data(:,t_els,it);
         dtmp_cfy = -sdata(6).data(:,t_els,it);
      
         [xsort ind] = sort(dtmpx(:));
         ysort = dtmpy(ind);
         cfx = dtmp_cfx(ind);
         cfy = dtmp_cfy(ind);    

         % Rotate imported values according to simulation time   
         t_sim = tstamps(it);
         dtheta = ptch_amp*pi/180*sin(omega*(t_sim-ptch_start)+phase)-ptch_amp*pi/180*sin(phase);
         theta = atan2(sty_top,stx_top);
         theta_new = theta-dtheta;        % emperically decided sign
         sty_new = sin(theta_new);
         stx_new = cos(theta_new);

         xnew = xt_imp;
         ynew = yt_imp;

         % positive clockwise              
         rot = [cos(dtheta) sin(dtheta); ...
                -sin(dtheta) cos(dtheta)];

         coords = rot*[transpose(xnew)-axis_x0; transpose(ynew)-axis_y0];
         xrnew = coords(1,:) + axis_x0;
         yrnew = coords(2,:) + axis_y0;
         %% end of rotation   

                         
         x_min = min(xsort);
         x_max = max(xsort);

         [val ind1] = min(abs(xrnew-x_min));
         x_ref = xrnew(ind1:ind1+length(xsort)-1);
         y_ref = yrnew(ind1:ind1+length(xsort)-1);       

         stx_ref = stx_new(ind1:ind1+length(xsort)-1);
         sty_ref = sty_new(ind1:ind1+length(xsort)-1);       

         cf = cfx.*(stx_ref') + cfy.*(sty_ref');

%        Find transition location
         ind2 = find(xsort>0.1);
         trxsort = xsort(ind2);
         trcf = cf(ind2);            

         order=1;
         framelen=51;      
         trcf1 = sgolayfilt(trcf,order,framelen);
         cfgrad = gradient(trcf,trxsort);

         order=1;
         framelen=21;      
         cfgrad1 = sgolayfilt(cfgrad,order,framelen);
         trind = find(abs(cfgrad1)>0.035,1);
         trloc = trxsort(trind);

         xtr = [xtr trloc];   

         if (ifplot)   
           pvar = plot(xsort,cf, 'b.', 'MarkerSize', 6);
%           pvar = plot(trxsort,cfgrad1, 'b.', 'MarkerSize', 6);
           grid on
           hold on 
%           xlim([0.05 .15])    
%           set(gca,'Ydir', 'reverse')
           ylim([-0.0035 0.005]);

           ylims = get(gca,'YLim');
           trplot = plot([trloc trloc], ylims, '--r'); 
%           grid on   
%           hold on
           lgs{1} =  ['T=' num2str(tstamps(it))]; 
           lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'SouthWest', 'Fontsize', lfs, 'Box', 'off');
%           if nplots == 0 
             ylabel('C_{f}', 'Interpreter', 'tex', 'Fontsize', fs);
             xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
%           end
           hold off 

%          surf_x = [surf_x; dtmpx(:)'/Chord];
%          surf_t = [surf_t; tstamps(it)*ones(1,length(dtmpx(:)))];
%          surf_v = [surf_v; dtmp_v(:)'];
         end    

         surf_x = [surf_x; xsort'/Chord];
         surf_t = [surf_t; tstamps(it)*ones(1,length(xsort))];
         surf_v = [surf_v; cf'];
         surf_c = [surf_c; sign(cf)'];

         ind2 = find(cf<0,1);
         xbub_st = xsort(ind2);


         nplots = nplots+1;   
%         mov(nplots) = getframe(gcf);

         svfname = sprintf('%0.5d', nplots);   
         svfname = ['pitch750k' svfname '.png'];
         destn = 'plots/';   
%         SaveFig(gcf, svfname, destn, 1)
      end
      tlast = tstamps(it);
      if (ifplot)
        pause(0.01)
      end
    end
  end
end           

figure(2)
%ax1 = subplot(1,2,1);
ax1=axes;
splot=surf(ax1,surf_x,(surf_t-ptch_start)/Tosc-0.0,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
colorbar;
view(2)
ylabel('$\frac{t}{T_{osc}}$', 'Interpreter','Latex', 'rot', 0, 'FontSize', 16)
xlabel('x/c', 'FontSize', 16)
xlim([0 1])
%ylim([min(surf_t(:)) max(surf_t(:))])
axis tight
hold on
set(ax1, 'YTickMode', 'manual')
yticks = [1 1.25 1.5 1.75 2 2.25 2.50]-1.00;
set(ax1, 'YTick', yticks);
set(ax1, 'FontSize', 14)
%set(ax1, 'PlotBoxAspectRatio', [1 1.5 1])

svfname = ['cf_time_surf750k.eps'];
destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)



%cplot=contour(ax1,surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
%view(2)
%svfname = ['cf_time_surf_contour.eps'];
%destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)

figure(3)
%ax2=subplot(1,2,2);
ax3=axes;
gplot=surf(ax3,surf_x,(surf_t-ptch_start)/Tosc-0.0,surf_c,surf_c,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
set(ax3,'Color', 'none')
view(2)
colormap(ax3,'gray');
%ylabel('t/T_{osc}', 'FontSize', 16)
xlabel('x/c', 'FontSize', 16)
xlim([0 1])
axis tight
set(ax3, 'YTickMode', 'manual')
yticks = [0:7]*0.25;
set(ax3, 'YTick', yticks);
set(ax3, 'FontSize', 14)

ax4=axes;
gplot2=surf(ax4,surf_x,(surf_t-ptch_start)/Tosc,surf_c,surf_c,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); hold on
view(2)
axis tight
lncol1 = 'blue';
lncol2 = 'red';
ln1 = line([0 1], [1.25 1.25]-1.0, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', lncol1, 'Parent', ax4);
ln2 = line([0 1], [1.50 1.50]-1.0, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', lncol2, 'Parent', ax4);
ln3 = line([0 1], [1.75 1.75]-1.0, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', lncol1, 'Parent', ax4);
ln4 = line([0 1], [2.00 2.00]-1.0, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', lncol2, 'Parent', ax4);
ln4 = line([0 1], [2.25 2.25]-1.0, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', lncol1, 'Parent', ax4);
ln4 = line([0 1], [2.50 2.50]-1.0, [2 2], 'LineStyle', '--', 'LineWidth', 1.0, 'Color', lncol2, 'Parent', ax4);

set(ax4,'YAxisLocation', 'right');
set(ax4, 'XTick', []);
set(ax4, 'Color', 'none')
set(ax4, 'YTickMode', 'manual')
yticks = [0:7]*0.25;
set(ax4, 'YTick', yticks);
ylbl = {'-p/2', '0', 'p/2', 'p', '3p/2', '2p', '5p/2'};
set(ax4, 'YTickLabel', ylbl, 'FontName', 'symbol', 'FontSize', 14);
%ylabel('Oscillation phase', 'Font','Helvetica', 'FontSize', 16)
%set(ax4, 'PlotBoxAspectRatio', [1 1.5 1])

%colorbar
set(ax3,'OuterPosition', get(ax1,'OuterPosition'));
set(ax4,'OuterPosition', get(ax1,'OuterPosition'));

set(ax3,'Position', get(ax1,'Position'));
set(ax4,'Position', get(ax1,'Position'));

figure(3)
svfname = ['cf_time_surf_grey750k.eps'];
destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)


figure(4)
ax5=axes;
pitch_time = surf_t(:,1);
alpha0 = 3.4;
pitch_alpha = alpha0 + ptch_amp*sin(omega*(pitch_time-ptch_start)+phase);
trplot=plot(ax5,pitch_alpha,xtr, '.'); hold on
xlabel('\alpha^{o}')
ylabel('x/c')
title('Transition phase plot')
svfname = ['transition_phase750k.eps'];
destn = 'plots/';   
SaveFig(gcf, svfname, destn, 1)


% ncontours = 2;
% cont_vec = linspace(min(surf_v(:)),0,ncontours);
% cplot=contour(surf_x,surf_t,surf_v,cont_vec);

%mov2 = mov(1:nplots); 
%movie2avi(mov2, 'cp_movie.avi', 'compression', 'None', 'fps', 15, 'quality', 75);
%      for ies=1:selt
%      %     scatter(x(:,i),sdata(1).data(:,i,1), '.')
%      %     scatter(x(:,i),sdata(1).data(:,i,maxtsaves), 'd')
%           scatter(sdata(1).data(:,ies,maxtsaves),sdata(3).data(:,ies,maxtsaves), 'd')
%      
%      end
%      set(gca,'Ydir','reverse')
     



