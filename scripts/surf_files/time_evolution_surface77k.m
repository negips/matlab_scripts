% Plot cp files

clear
clc
close all

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 'naca0012_re77k/';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
lafs = 20;              % latex fontsize
ifcols = 1;
ifplot = 0;
tlast = 19.05;
tstart0 = tlast; 
destn = 'plots/';
ifcontour=0;
iftr=0;                 % Plot transition location contour
iftrabs=0;              % plot point of absolute instability
iftrgrey=0;             % Plot transition location contour
ifsave = 0;             % if save figures
ifdatasave=1;           % if save mat file
axfs=15;                % axis font size

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
maxframes = nfiles*100;

if (ifplot)
  h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);
end

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

% Load mean aoa normals for polynomial order 9
%snormals = importdata('surf_normals.67.N5');
%
%x_imp5 = snormals.data(:,1);
%[x_imp5 I] = sort(x_imp5);
%y_imp5 = snormals.data(I,2);
%
%snx5 = -snormals.data(I,4);
%sny5 = -snormals.data(I,5);
%
%ind = sny5>0;
%snx_top5 = snx5(find(ind));
%sny_top5 = sny5(find(ind));
%xt_imp5  = x_imp5(find(ind));
%yt_imp5  = y_imp5(find(ind));
%
%snx_bot5 = snx5(find(~ind));
%sny_bot5 = sny5(find(~ind));
%xb_imp5  = x_imp5(find(~ind));
%yb_imp5  = y_imp5(find(~ind));
%
%for i=1:length(snx_top5)
%  stx_top5(i) = sny_top5(i);
%  sty_top5(i) = -snx_top5(i);
%end
%
%for i=1:length(snx_bot5)
%  stx_bot5(i) = -sny_bot5(i);
%  sty_bot5(i) = snx_bot5(i);
%end
%---------------------------------------- 
snormals = importdata('surf_normals.0012.N9');

x_imp9 = snormals.data(:,1);
[x_imp9 I] = sort(x_imp9);
y_imp9 = snormals.data(I,2);

snx9 = -snormals.data(I,4);
sny9 = -snormals.data(I,5);

ind = sny9>0;
snx_top9 = snx9(find(ind));
sny_top9 = sny9(find(ind));
xt_imp9  = x_imp9(find(ind));
yt_imp9  = y_imp9(find(ind));

snx_bot9 = snx9(find(~ind));
sny_bot9 = sny9(find(~ind));
xb_imp9  = x_imp9(find(~ind));
yb_imp9  = y_imp9(find(~ind));

for i=1:length(snx_top9)
  stx_top9(i) = sny_top9(i);
  sty_top9(i) = -snx_top9(i);
end

for i=1:length(snx_bot9)
  stx_bot9(i) = -sny_bot9(i);
  sty_bot9(i) = snx_bot9(i);
end
%% 

U0=1.;
kred=0.5;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=2*pi/omega;
ptch_amp = 1.3;
ptch_start = 0.;
axis_x0 = 0.18;
axis_y0 = 0.0;
phase=0;

nplots = 0;
bcnt = 0;
surf_x5 = [];
surf_t5 = [];
surf_v5 = [];
surf_c5 = [];
surf_p5 = [];
npts5 = 0;

surf_x9 = [];
surf_t9 = [];
surf_v9 = [];
surf_c9 = [];
surf_p9 = [];
npts9 = 0;

fsi_io1 = importdata('fsi_io.108');
fsi_io2 = importdata('fsi_io.120');

icalld = 0;
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    if (lx1(1)==6)
      sty_top = sty_top5;
      stx_top = stx_top5;
      xt_imp = xt_imp5;
      yt_imp = yt_imp5;
    elseif (lx1(1)==10)    
      sty_top = sty_top9;
      stx_top = stx_top9;
      xt_imp = xt_imp9;
      yt_imp = yt_imp9;
    end  

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

      for i=1:2   % 2 leading elements
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
         end   
         dtmpx = sdata(1).data(:,t_els,it);
         dtmpy = sdata(2).data(:,t_els,it);
%         dtmp_v = sdata(3).data(:,t_els,it);    
         dtmp_v = -sdata(5).data(:,t_els,it);
         dtmp_cfx = -sdata(5).data(:,t_els,it);
         dtmp_cfy = -sdata(6).data(:,t_els,it);
         dtmp_cp = sdata(3).data(:,t_els,it); 
      
         [xsort ind] = sort(dtmpx(:));
         ysort = dtmpy(ind);
         cfx = dtmp_cfx(ind);
         cfy = dtmp_cfy(ind);
         cp = dtmp_cp(ind);    

         % Rotate imported values according to simulation time   
         t_sim = tstamps(it);
%         dtheta = ptch_amp*pi/180*sin(omega*(t_sim-ptch_start)+phase)-ptch_amp*pi/180*sin(phase);
         if (t_sim<108.00028)   
           dtheta = interp1(fsi_io1.data(:,2),fsi_io1.data(:,4),t_sim);
         else
           dtheta = interp1(fsi_io2.data(:,2),fsi_io2.data(:,4),t_sim);              
         end   
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

         if (ifplot)   
           pvar = plot(xsort,cf, 'b.', 'MarkerSize', 6);
           grid on
%           xlim([0.05 .30])    
%           set(gca,'Ydir', 'reverse')
           ylim([-0.005 0.015]);
%           grid on   
%           hold on
           lgs{1} =  ['T=' num2str(tstamps(it))]; 
           lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'SouthWest', 'Fontsize', lfs, 'Box', 'off');
%           if nplots == 0 
             ylabel('$\tau_{w}$', 'Interpreter', 'latex', 'Fontsize', fs);
             xlabel('$x/C$', 'Interpreter', 'latex', 'Fontsize', fs);
%           end   

%          surf_x = [surf_x; dtmpx(:)'/Chord];
%          surf_t = [surf_t; tstamps(it)*ones(1,length(dtmpx(:)))];
%          surf_v = [surf_v; dtmp_v(:)'];
         end    

         if (lx1(1)==6)   
           surf_x5 = [surf_x5; xsort'/Chord];
           surf_t5 = [surf_t5; tstamps(it)*ones(1,length(xsort))];
           surf_v5 = [surf_v5; cf'];
           surf_c5 = [surf_c5; sign(cf)'];
           surf_p5 = [surf_p5; cp'];
           npts5=npts5+1;
         elseif(lx1(1)==10)    
           surf_x9 = [surf_x9; xsort'/Chord];
           surf_t9 = [surf_t9; tstamps(it)*ones(1,length(xsort))];
           surf_v9 = [surf_v9; cf'];
           surf_c9 = [surf_c9; sign(cf)'];
           surf_p9 = [surf_p9; cp'];
           npts9=npts9+1;
         end   

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
end  % end of files 


if ifdatasave
  save('re77k_surface.mat')
end

%splot_100k

     



