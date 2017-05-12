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
tlast = 19.0;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

nplots = 0;
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
      xtels_mean = mean(sdata(1).data(:,t_els,1));
      [val ind] = min(xtels_mean);
      t_els(ind)=[];                  % remove leading edge element 

      xtels_mean = mean(sdata(1).data(:,t_els,1));
      [val ind] = min(xtels_mean);
      t_els(ind)=[];                  % remove 2nd leading edge element 
    
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
         [xsort ind] = sort(dtmpx(:));
          ysort = dtmpy(ind);
          vsort = dtmp_v(ind);  

         surf_x = [surf_x; xsort'/Chord];
         surf_t = [surf_t; tstamps(it)*ones(1,length(xsort(:)))];
         surf_v = [surf_v; vsort'];      
         surf_c = [surf_c; sign(vsort)'];      

%         figure(h1)  
%         pvar = plot(xsort/Chord, vsort, 'b.', 'MarkerSize', 8);
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
% ncontours = 2;
% cont_vec = linspace(min(surf_v(:)),0,ncontours);
% cplot=contour(surf_x,surf_t,surf_v,cont_vec);
svfname = ['cf_time_surf.eps'];
destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)

figure(3)
splot=surf(surf_x,surf_t,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.25); hold on
splot=surf(surf_x,surf_t,surf_c,surf_c,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp');
%colorbar;
view(2)
ylabel('Time')
xlabel('x')
xlim([0 1])
colormap('gray')
svfname = ['cf_time_surf_grey.eps'];
destn = 'plots/';   
%SaveFig(gcf, svfname, destn, 1)


%mov2 = mov(1:nplots); 
%movie2avi(mov2, 'cp_movie.avi', 'compression', 'None', 'fps', 15, 'quality', 75);
%      for ies=1:selt
%      %     scatter(x(:,i),sdata(1).data(:,i,1), '.')
%      %     scatter(x(:,i),sdata(1).data(:,i,maxtsaves), 'd')
%           scatter(sdata(1).data(:,ies,maxtsaves),sdata(3).data(:,ies,maxtsaves), 'd')
%      
%      end
%      set(gca,'Ydir','reverse')
     



