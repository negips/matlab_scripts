% Plot cp files

clear
clc
close all

% addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 're750k_pitch';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
lafs = 24;              % Latex font size
ifcols = 1;
destn = 'plots/';

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

%% Load mean aoa normals for polynomial order 5
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
%% 




[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 31.0;
tmax = 39.00001;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

cfavg = [];
cfavgx = [];
cfavgy = [];
ncf_pts = 0;
cf_start=41.25;
cf_end = 48.04;
ifcfplot = 0;
nplots = 0;

cols=lines(50);

for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};

    ndim=3;    
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr,ndim);

    if (tstamps(1)>tmax)
       break
    end     

    if (lx1(1)==6)
      sty_top = sty_top5;
      stx_top = stx_top5;
      x_imp = x_imp5;
      y_imp = y_imp5;
    elseif (lx1(1)==7)    
      sty_top = sty_top6;
      stx_top = stx_top6;
      x_imp = x_imp6;
      y_imp = y_imp6;
    elseif (lx1(1)==9)    
      sty_top = sty_top8;
      stx_top = stx_top8;
      x_imp = x_imp8;
      y_imp = y_imp8;
    end  

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

    xmin=min(x(:));
    xmax=max(x(:));
    Chord=xmax-xmin;
%    Chord=1; 

    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast)
         
         if (tstamps(it)>tmax)
            break
         end         

         if nplots>0
           delete(pvar)
           delete(pvar2)
         end   
         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);
         dtmp_v = sdata(3).data(:,:,it);    
%         dtmp_v = -sdata(5).data(:,:,it);    
         cf = -sdata(5).data(:,:,it);
         cp = sdata(3).data(:,:,it);

         [xsort ind] = sort(dtmpx(:));
         ysort = dtmpy(:);
         ysort = ysort(ind);
         cf = cf(:);
         cf = cf(ind);
         cp = cp(ind); 

%%        Average cf   
         if tstamps(it)>cf_start && tstamps(it)<cf_end

           ncf_pts = ncf_pts+1;
           avgbeta = 1/ncf_pts;
           avgalpha = 1 - avgbeta;

           if avgalpha==0
             cfavg = cf;
             cfavgx = xsort/Chord;
             cfavgy = ysort/Chord;
           else               
             cfavg =  avgalpha*cfavg +  avgbeta*cf;
             cfavgx =  avgalpha*cfavgx +  avgbeta*xsort/Chord;
             cfavgy =  avgalpha*cfavgy +  avgbeta*ysort/Chord;
           end

         end

         if (tstamps(it)>cf_end) && (ifcfplot==0)
           cfplot = plot(cfavgx,cfavg, ' .r');
           ifcfplot=1; 
         end    

%%
         % Rotate imported values according to simulation time   
         t_sim = tstamps(it);
         dtheta = ptch_amp*pi/180*sin(omega*(t_sim-ptch_start)+phase)-ptch_amp*pi/180*sin(phase);
         alphamin = 2.4*pi/180;
         alpha = alphamin + dtheta;
%         theta = atan2(sty_top,stx_top);
%         theta_new = theta-dtheta;        % emperically decided sign
%         sty_new = sin(theta_new);
%         stx_new = cos(theta_new);

         xnew = x_imp;
         ynew = y_imp;

         % positive clockwise              
         rot = [cos(dtheta) sin(dtheta); ...
                -sin(dtheta) cos(dtheta)];

         coords = rot*[transpose(xnew)-axis_x0; transpose(ynew)-axis_y0];
         xrnew = coords(1,:) + axis_x0;
         yrnew = coords(2,:) + axis_y0;
         %% end of rotation   

%%
         
         figure(h1)      
%         pvar = plot(x(:)/Chord,dtmp_v(:), 'b.', 'MarkerSize', 10);
%         pvar = plot(xsort/Chord,cp, ' .', 'Color', cols(it,:), 'MarkerSize', 15);
         pvar = plot(xsort,ysort, 'b.', 'MarkerSize', 15); hold on
         pvar2 = plot(xrnew,yrnew, '.r', 'MarkerSize', 15);
        
%         set(gca,'Ydir', 'reverse')
%         ylim([-3.5 1.1]);
%         xlim([-0.01 1.000])
%         ylim([-0.004 0.0060])
         xlim([0.96 1])
         ylim([-0.1 -0.05])
         grid on   
         hold on
         lgs{1} =  ['T=' num2str(tstamps(it)) '; alpha=' num2str(round(100*alpha*180/pi)/100)]; 
         lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'North', 'Fontsize', lfs, 'Box', 'off');
         if nplots == 0 
           ylabel('$\tau_{w}$', 'Interpreter', 'latex', 'Fontsize', lafs);
           xlabel('$x/C$', 'Interpreter', 'latex', 'Fontsize', lafs);
         end   
         nplots = nplots+1;   
%         mov(nplots) = getframe(gcf);

         svfname = sprintf('%0.4d', nplots);   
         svfname = ['cf_fig' svfname '.eps'];
         destn = 'plots/';   
%         SaveFig(gcf, svfname, destn, 1)
      end
      tlast = tstamps(it);
      pause(0.01)
    end           % it
  end
end           

%mov2 = mov(1:nplots); 
%movie2avi(mov2, 'cp_movie.avi', 'compression', 'None', 'fps', 15, 'quality', 75);
%      for ies=1:selt
%      %     scatter(x(:,i),sdata(1).data(:,i,1), '.')
%      %     scatter(x(:,i),sdata(1).data(:,i,maxtsaves), 'd')
%           scatter(sdata(1).data(:,ies,maxtsaves),sdata(3).data(:,ies,maxtsaves), 'd')
%      
%      end
%      set(gca,'Ydir','reverse')
     



