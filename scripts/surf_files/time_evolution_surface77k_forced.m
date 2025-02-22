% Plot cp files

clear
clc
close all

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

U0=1.;
kred=0.0477;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=2*pi/omega;
ptch_amp = 5.5;
ptch_start = 16.0;
axis_x0 = 0.186;
axis_y0 = 0.0;
phase=0;


fol = 'naca0012_re77k_k0.0477/';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
lafs = 20;              % latex fontsize
ifcols = 1;
ifplot = 0;
tlast = 16.00;
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
snormals = importdata('surf_normals.0012.N7');
x_imp7 = snormals.data(:,1);
[x_imp7 I] = sort(x_imp7);
y_imp7 = snormals.data(I,2);

snx7 = -snormals.data(I,4);
sny7 = -snormals.data(I,5);

ind = sny7>0;
snx_top7 = snx7(find(ind));
sny_top7 = sny7(find(ind));
xt_imp7  = x_imp7(find(ind));
yt_imp7  = y_imp7(find(ind));

snx_bot7 = snx7(find(~ind));
sny_bot7 = sny7(find(~ind));
xb_imp7  = x_imp7(find(~ind));
yb_imp7  = y_imp7(find(~ind));

for i=1:length(snx_top7)
  stx_top7(i) = sny_top7(i);
  sty_top7(i) = -snx_top7(i);
end

for i=1:length(snx_bot7)
  stx_bot7(i) = -sny_bot7(i);
  sty_bot7(i) = snx_bot7(i);
end
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

nplots = 0;
bcnt = 0;
surf_x7 = [];
surf_t7 = [];
surf_v7 = [];
surf_c7 = [];
surf_p7 = [];
npts7 = 0;

surf_x9 = [];
surf_t9 = [];
surf_v9 = [];
surf_c9 = [];
surf_p9 = [];
npts9 = 0;

%fsi_io1 = importdata('fsi_io.108');
%fsi_io2 = importdata('fsi_io.120');

icalld = 0;
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr,ndim);

    if (lx1(1)==8)
      sty_top = sty_top7;
      stx_top = stx_top7;
      xt_imp = xt_imp7;
      yt_imp = yt_imp7;
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
%           delete(pvar2)
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
         dtheta = ptch_amp*pi/180*sin(omega*(t_sim-ptch_start)+phase);

         % positive clockwise              
         Rot = [cos(dtheta) sin(dtheta); ...
               -sin(dtheta) cos(dtheta)]; 

%         Rot2 = [cos(dtheta) -sin(dtheta); ...
%                sin(dtheta) cos(dtheta)]; 
              
         st_rot = Rot*[stx_top; sty_top];   
         stx_new = st_rot(1,:);
         sty_new = st_rot(2,:);

%         st_rot2 = Rot2*[stx_top; sty_top];   
%         stx_new2 = st_rot2(1,:);
%         sty_new2 = st_rot2(2,:);

         xnew = xt_imp;
         ynew = yt_imp;

         coords = Rot*[transpose(xnew)-axis_x0; transpose(ynew)-axis_y0];
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

%         stx_ref2 = stx_new2(ind1:ind1+length(xsort)-1);
%         sty_ref2 = sty_new2(ind1:ind1+length(xsort)-1);       
%
%         cf2 = cfx.*(stx_ref2') + cfy.*(sty_ref2');

         if (ifplot)   
           pvar = plot(xsort,cf, 'b.', 'MarkerSize', 6);
           grid on
%           pvar2 = plot(xsort,cf2, 'r.', 'MarkerSize', 6); hold off

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

         if (lx1(1)==8)   
           surf_x7 = [surf_x7; xsort'/Chord];
           surf_t7 = [surf_t7; tstamps(it)*ones(1,length(xsort))];
           surf_v7 = [surf_v7; cf'];
           surf_c7 = [surf_c7; sign(cf)'];
           surf_p7 = [surf_p7; cp'];
           npts7=npts7+1;
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
%         svfname = ['pitch750k' svfname '.png'];
%         destn = 'plots/';   
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
  save(['re77k_k' num2str(kred*10000) '_surface.mat'])
end

%splot_100k

     



