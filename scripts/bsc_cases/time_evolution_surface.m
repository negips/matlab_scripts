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
tlast = 18.85;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

% Load mean aoa normals
snormals = importdata('surf_normals.67');

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
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    xmax = max(x(:));
    xmin = min(x(:));      
    Chord=xmax-xmin;

%   trying to find top and bottom elements        
    [val ind] = max(sdata(1).data(:,:,1));
    n12=find(ind==12);
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
         if nplots>0
%           delete(pvar)
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
         dtheta = ptch_amp*pi/180*sin(omega*(t_sim-ptch_start));
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

%         pvar = plot(xsort,cf, 'b.', 'MarkerSize', 6);
%         xlim([0.05 .15])    
%         set(gca,'Ydir', 'reverse')
%         ylim([-3.5 1.1]);
%         grid on   
%         hold on
%         lgs{1} =  ['T=' num2str(tstamps(it))]; 
%         lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'North', 'Fontsize', lfs, 'Box', 'off');
%         if nplots == 0 
%           ylabel('C_{p}', 'Interpreter', 'tex', 'Fontsize', fs);
%           xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
%         end   

%         surf_x = [surf_x; dtmpx(:)'/Chord];
%         surf_t = [surf_t; tstamps(it)*ones(1,length(dtmpx(:)))];
%         surf_v = [surf_v; dtmp_v(:)'];

         surf_x = [surf_x; xsort'/Chord];
         surf_t = [surf_t; tstamps(it)*ones(1,length(xsort))];
         surf_v = [surf_v; cf'];
         surf_c = [surf_c; sign(cf)'];

         ind2 = find(cf<0,1);
         xbub_st = xsort(ind2);

         if xbub_st<0.18                  % arbitrarily set aposteriory
           bcnt = bcnt+1; 
           bubble_start(bcnt) = xbub_st; 

           tmp_x = xsort(ind2:end);
           tmp_cf = cf(ind2:end);
           ind3 = find(tmp_cf>0,1);           % first time cf goes positive
           bubble_end(bcnt) = tmp_x(ind3);
           bubble_time(bcnt) = tstamps(it);
         end    

         nplots = nplots+1;   
%         mov(nplots) = getframe(gcf);

         svfname = sprintf('%0.4d', nplots);   
         svfname = ['cp_fig' svfname '.png'];
         destn = 'plots/';   
%         SaveFig(gcf, svfname, destn, 1)
      end
      tlast = tstamps(it);
      pause(0.01)
    end
  end
end           

figure(2)
splot=surf(surf_x,surf_t,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
colorbar;
view(2)
ylabel('Time')
xlabel('x')
xlim([0 1])
hold on
svfname = ['cf_time_surf.eps'];
destn = 'plots/';   
SaveFig(gcf, svfname, destn, 1)

cplot=contour(surf_x,surf_t,surf_v,[0 0], 'LineColor', 'k', 'LineWidth', 1.5);
svfname = ['cf_time_surf_contour.eps'];
destn = 'plots/';   
SaveFig(gcf, svfname, destn, 1)

figure(3)
ax2=axes;
set(ax2,'Color', 'none')
gplot=surf(ax2,surf_x,surf_t,surf_c,surf_c,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
view(2)
colormap(ax2,'gray');
ylabel('Time')
xlabel('x')
xlim([0 1])

svfname = ['cf_time_surf_grey.eps'];
destn = 'plots/';   
SaveFig(gcf, svfname, destn, 1)


figure(4)
plot(bubble_time,bubble_start,'b'); hold on
plot(bubble_time,bubble_end, 'r')


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
     



