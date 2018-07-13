% Plot cp files

clear
clc
close all

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 're100k';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 18.0;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition', [0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

% Load mean aoa normals
[xgll wts] = lglnodes(5);

snormals = importdata('surf_normals.67.N5');

x_imp = snormals.data(:,1);
%[x_imp I] = sort(x_imp);
I = [1:length(x_imp)];
y_imp = snormals.data(I,2);

snx = snormals.data(I,4);
sny = snormals.data(I,5);

ind = sny<0;
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

x_imp2 = reshape(x_imp,6,length(x_imp)/6);
y_imp2 = reshape(y_imp,6,length(y_imp)/6);
snx_imp = reshape(snx,6,length(x_imp)/6);
sny_imp = reshape(sny,6,length(y_imp)/6);

%% 
% Load transition point calculations

load('tr.mat')

%%
U0=1.;
kred=0.5;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
ptch_amp = 1.3;
ptch_start = 0.;
axis_x0 = 0.35;
axis_y0 = 0.034;

nplots = 0;
bcnt = 0;
surf_x = [];
surf_t = [];
surf_v = [];
surf_c = [];
icalld = 0;
cl_hist = [];
time_hist = [];
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    if timeout>19.0
      break
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

      for i=1:2   % 5 leading elements
        xtels_mean = mean(sdata(1).data(:,t_els,1));
        [val ind] = min(xtels_mean);
        t_els(ind)=[];                  % remove leading edge element

        xbels_mean = mean(sdata(1).data(:,b_els,1));
        [val ind] = min(xbels_mean);
        b_els(ind)=[];                  % remove leading edge element
      end 
    
      total_top_els = length(t_els);
      total_bot_els = length(b_els);
      icalld = icalld+1;
    else
      ntop = length(t_els);
      while ntop>total_top_els
        xtels_mean = mean(sdata(1).data(:,t_els,1));
        [val ind] = min(xtels_mean);
        t_els(ind)=[];                  % remove leading edge element 
        ntop = length(t_els);
      end

      nbot = length(b_els);
      while nbot>total_bot_els
        xbels_mean = mean(sdata(1).data(:,b_els,1));
        [val ind] = min(xbels_mean);
        b_els(ind)=[];                  % remove leading edge element 
        nbot = length(t_els);
      end


    end    

      
    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast)
         if nplots>0
%           delete(pvar)
         end   
         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);

         dtmp_v = sdata(3).data(:,:,it);
         cl_p = sdata(3).data(:,:,it);
         cl_v = sdata(6).data(:,:,it);      

%         dtmp_v = -sdata(5).data(:,t_els,it);
%         dtmp_cfx = sdata(5).data(:,t_els,it);
%         dtmp_cfy = sdata(6).data(:,t_els,it);

         % Rotate imported values according to simulation time   
         t_sim = tstamps(it);
         dtheta = ptch_amp*pi/180*sin(omega*(t_sim-ptch_start));
         theta = atan2(sty_top,stx_top);
         theta_new = theta-dtheta;        % emperically decided sign
         sty_new = sin(theta_new);
         stx_new = cos(theta_new);

         xnew = x_imp;
         ynew = y_imp;

         % positive clockwise              
         rot = [cos(dtheta) sin(dtheta); ...
                -sin(dtheta) cos(dtheta)];

         coords = rot*[transpose(xnew)-axis_x0; transpose(ynew)-axis_y0];
         xrnew = coords(1,:) + axis_x0;
         yrnew = coords(2,:) + axis_y0;
         %% end of rotation

         xrnew=reshape(xrnew,6,length(xrnew)/6);   
         yrnew=reshape(yrnew,6,length(yrnew)/6);   

         snx_new = zeros(lx1(1),selt);
         sny_new = zeros(lx1(1),selt);
         indices = [];
         resid = [];
         cl = 0;    
%         for jj=1:selt
%           for kk=1:selt
%              if ( (abs(dtmpx(1,jj)-xrnew(1,kk))<1e-6) &&  (abs(dtmpx(lx1(1),jj)-xrnew(lx1(1),kk))<1e-6) )
%                indices = [indices kk];
%                resid = [resid  sqrt((dtmpx(1,jj)-xrnew(1,kk))^2 + (dtmpx(lx1(1),jj)-xrnew(lx1(1),kk))^2)];
%                snx_new(:,jj)=snx_imp(:,kk); 
%                sny_new(:,jj)=sny_imp(:,kk);
%                cl_el = (cl_p(:,jj) + 0*cl_v(:,jj)).*sny_new(:,jj);
%                intg = cl_el.*wts*abs(( max(dtmpx(:,jj)) - min(dtmpx(:,jj)) ))/2;  
%                cl=cl+sum(intg); 
%                break
%              end
%           end
%         end
%         cl_hist = [cl_hist cl];
%         time_hist = [time_hist tstamps(it)]; 


%         pvar = plot(xsort,cf, 'b.', 'MarkerSize', 6);
%         xlim([0.05 .15])    
%         set(gca,'Ydir', 'reverse')
%         ylim([-3.5 1.1]);
%         grid on   
%         hold on
%         tr_pt_uv = interp1(tr_time,trx_uv,tstamps(it));
%         tr_pt_ww = interp1(tr_time,trx_ww,tstamps(it));
%         gca_ylims = get(gca,'Ylim');   
%         ptruv = plot([tr_pt_uv tr_pt_uv], gca_ylims, '--b'); 
%         ptrww = plot([tr_pt_ww tr_pt_ww], gca_ylims, '--r');   

%         surf_x = [surf_x; xsort'/Chord];
%         surf_t = [surf_t; tstamps(it)*ones(1,length(xsort))];
%         surf_v = [surf_v; cf'];
%         surf_c = [surf_c; sign(cf)'];

         nplots = nplots+1;   
%         mov(nplots) = getframe(gcf);

         svfname = sprintf('%0.4d', nplots);   
         svfname = ['cp_fig' svfname '.png'];
         destn = 'plots/';   
%         SaveFig(gcf, svfname, destn, 1)

      end                                 % if tstamps(it)>tlast
      tlast = tstamps(it);
      pause(0.01)
    end           % it
%    plot(tstamps,(sintegrals(:,2,3) + sintegrals(:,2,6))*6, 'r'); hold on
    cl_hist = [cl_hist ;(sintegrals(:,2,3)+sintegrals(:,2,6))];
    time_hist = [time_hist; tstamps];
  end
end           

fac=7.6621
figure(1)
plot(time_hist, cl_hist*fac)

%figure(2)
%splot=surf(surf_x,surf_t,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
%colorbar;
%view(2)
%ylabel('Time')
%xlabel('x')
%xlim([0 1])
%hold on
%svfname = ['cf_time_surf.eps'];
%destn = 'plots/';   
%%SaveFig(gcf, svfname, destn, 1)
%
%cplot=contour(surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
%svfname = ['cf_time_surf_contour.eps'];
%destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)



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
     



