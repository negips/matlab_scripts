% Plot cp files

clear
clc
close all

%addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 're750k_pitch';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
ifplot = 0;             % plot individual wall profiles
tlast = 46.05;           % start from this time
tend = 100;             % stop at this time
destn = 'plots/';
ifcp = 0;               % plot pressure instead of cf
ifdatasave=1;           % save data into a mat file
  datafile='re750k_surfaceN9.mat';
lafs = 22;              % Latex font size

U0=1.;
kred=0.4;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=2*pi/omega;
%Tosc=1;
ptch_amp = 1.0;
ptch_start = 6.;
axis_x0 = 0.35;
axis_y0 = 0.034;
phase=-pi/2;

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
maxframes = nfiles*100;

if (ifplot)
  h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);
end

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

% Load mean aoa normals for polynomial order 5
snormals = importdata('surf_normals.24.N5');

x_imp5 = snormals.data(:,1);
[x_imp5 I] = sort(x_imp5);
y_imp5 = snormals.data(I,2);

snx5 = -snormals.data(I,4);
sny5 = -snormals.data(I,5);

ind = sny5>0;
snx_top5 = snx5(find(ind));
sny_top5 = sny5(find(ind));
xt_imp5  = x_imp5(find(ind));
yt_imp5  = y_imp5(find(ind));

snx_bot5 = snx5(find(~ind));
sny_bot5 = sny5(find(~ind));
xb_imp5  = x_imp5(find(~ind));
yb_imp5  = y_imp5(find(~ind));

for i=1:length(snx_top5)
  stx_top5(i) = sny_top5(i);
  sty_top5(i) = -snx_top5(i);
end

for i=1:length(snx_bot5)
  stx_bot5(i) = -sny_bot5(i);
  sty_bot5(i) = snx_bot5(i);
end
%---------------------------------------- 
snormals = importdata('surf_normals.24.N6');

x_imp6 = snormals.data(:,1);
[x_imp6 I] = sort(x_imp6);
y_imp6 = snormals.data(I,2);

snx6 = -snormals.data(I,4);
sny6 = -snormals.data(I,5);

ind = sny6>0;
snx_top6 = snx6(find(ind));
sny_top6 = sny6(find(ind));
xt_imp6  = x_imp6(find(ind));
yt_imp6  = y_imp6(find(ind));

snx_bot6 = snx6(find(~ind));
sny_bot6 = sny6(find(~ind));
xb_imp6  = x_imp6(find(~ind));
yb_imp6  = y_imp6(find(~ind));

for i=1:length(snx_top6)
  stx_top6(i) = sny_top6(i);
  sty_top6(i) = -snx_top6(i);
end

for i=1:length(snx_bot6)
  stx_bot6(i) = -sny_bot6(i);
  sty_bot6(i) = snx_bot6(i);
end
%---------------------------------------- 
snormals = importdata('surf_normals.24.N8');

x_imp8 = snormals.data(:,1);
[x_imp8 I] = sort(x_imp8);
y_imp8 = snormals.data(I,2);

snx8 = -snormals.data(I,4);
sny8 = -snormals.data(I,5);

ind = sny8>0;
snx_top8 = snx8(find(ind));
sny_top8 = sny8(find(ind));
xt_imp8  = x_imp8(find(ind));
yt_imp8  = y_imp8(find(ind));

snx_bot8 = snx8(find(~ind));
sny_bot8 = sny8(find(~ind));
xb_imp8  = x_imp8(find(~ind));
yb_imp8  = y_imp8(find(~ind));

for i=1:length(snx_top8)
  stx_top8(i) = sny_top8(i);
  sty_top8(i) = -snx_top8(i);
end

for i=1:length(snx_bot8)
  stx_bot8(i) = -sny_bot8(i);
  sty_bot8(i) = snx_bot8(i);
end
%---------------------------------------- 
snormals = importdata('surf_normals.24.N9');

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
surf_x5 = [];
surf_t5 = [];
surf_v5 = [];     % cf
surf_c5 = [];     % sign(cf)
surf_p5 = [];     % pressure

npts5 = 0;

surf_x8 = [];
surf_t8 = [];
surf_v8 = [];
surf_c8 = [];
surf_p8 = [];

npts8 = 0;

surf_x9 = [];
surf_t9 = [];
surf_v9 = [];
surf_c9 = [];
surf_p9 = [];

npts9 = 0;

cz = [];
cz_time = [];

icalld = 0;
for i = 1:nfiles
  if (tout(i)>=tlast) 
    fname = sfiles{i};
    
    ndim=3;
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr,ndim);

    if (lx1(1)==6)
      sty_top = sty_top5;
      stx_top = stx_top5;
      xt_imp = xt_imp5;
      yt_imp = yt_imp5;
    elseif (lx1(1)==7)    
      sty_top = sty_top6;
      stx_top = stx_top6;
      xt_imp = xt_imp6;
      yt_imp = yt_imp6;
    elseif (lx1(1)==9)    
      sty_top = sty_top8;
      stx_top = stx_top8;
      xt_imp = xt_imp8;
      yt_imp = yt_imp8;
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
      if (tstamps(it)>=tlast && tstamps(it)<=tend)
         if ifplot && nplots>0
           delete(pvar)
         end   
         dtmpx = sdata(1).data(:,t_els,it);
         dtmpy = sdata(2).data(:,t_els,it);
%         dtmp_v = sdata(3).data(:,t_els,it);
         dtmp_cp = sdata(3).data(:,t_els,it);
         pmax = max(dtmp_cp(:));

         dtmp_v = -sdata(5).data(:,t_els,it);
         dtmp_cfx = -sdata(5).data(:,t_els,it);
         dtmp_cfy = -sdata(6).data(:,t_els,it);
      
         [xsort ind] = sort(dtmpx(:));
         ysort = dtmpy(ind);
         cfx = dtmp_cfx(ind);
         cfy = dtmp_cfy(ind); 

         cp = dtmp_cp(ind);

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
         cpy = cp.*(sty_ref'); 

         if (ifplot)
           if (ifcp)
             pvar = plot(xsort,cp, 'b.', 'MarkerSize', 6);
             grid on
             set(gca,'Ydir', 'reverse')
             ylim([-1.1 0.1]);
             ylabel('C_{p}', 'Interpreter', 'tex', 'Fontsize', fs);
             xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
           else        
             pvar = plot(xsort,cf, 'b.', 'MarkerSize', 6);
             grid on
             ylim([-0.0035 0.005]);
%             if nplots == 0 
               ylabel('C_{f}', 'Interpreter', 'tex', 'Fontsize', fs);
               xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
%             end   
           end 
           lgs{1} =  ['T=' num2str(tstamps(it))]; 
           lg = legend(pvar,lgs, 'FontSize', lafs, 'Location', 'SouthWest', 'Box', 'off');

         end      % ifplot 

         if (lx1(1)==6)   
           surf_x5 = [surf_x5; xsort'/Chord];
           surf_t5 = [surf_t5; tstamps(it)*ones(1,length(xsort))];
           surf_v5 = [surf_v5; cf'];
           surf_c5 = [surf_c5; sign(cf)'];
           surf_p5 = [surf_p5; cp'];
           npts5=npts5+1;
         elseif(lx1(1)==7)    
           surf_x8 = [surf_x8; xsort'/Chord];
           surf_t8 = [surf_t8; tstamps(it)*ones(1,length(xsort))];
           surf_v8 = [surf_v8; cf'];
           surf_c8 = [surf_c8; sign(cf)'];
           surf_p8 = [surf_p8; cp'];
           npts8=npts8+1;
         elseif(lx1(1)==9)    
           surf_x8 = [surf_x8; xsort'/Chord];
           surf_t8 = [surf_t8; tstamps(it)*ones(1,length(xsort))];
           surf_v8 = [surf_v8; cf'];
           surf_c8 = [surf_c8; sign(cf)'];
           surf_p8 = [surf_p8; cp'];
           npts8=npts8+1;
         elseif(lx1(1)==10)    
           surf_x9 = [surf_x9; xsort'/Chord];
           surf_t9 = [surf_t9; tstamps(it)*ones(1,length(xsort))];
           surf_v9 = [surf_v9; cf'];
           surf_c9 = [surf_c9; sign(cf)'];
           surf_p9 = [surf_p9; cp'];
           npts9=npts9+1;
         end

         cz = [cz sintegrals(it,2,3)];
         cz_time = [cz_time tstamps(it)];   

         ind2 = find(cf<0,1);
         xbub_st = xsort(ind2);

%         if xbub_st<0.18                  % arbitrarily set aposteriory
%           bcnt = bcnt+1; 
%           bubble_start(bcnt) = xbub_st; 
%
%           tmp_x = xsort(ind2:end);
%           tmp_cf = cf(ind2:end);
%           ind3 = find(tmp_cf>0,1);           % first time cf goes positive
%           bubble_end(bcnt) = tmp_x(ind3);
%           bubble_time(bcnt) = tstamps(it);
%         end    

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


if ifdatasave
  disp(['Saving file ' datafile])
  save(datafile)
end

%splot_750k


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
     



