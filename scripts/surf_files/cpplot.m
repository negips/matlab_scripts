% Plot cp files

clear
clc
close all

% addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 're750k_aoa34';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
lafs = 24;              % Latex font size
ifcols = 1;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 2.70;
tmax = 500.00001;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

cfavg = [];
cfavgx = [];
cfavgy = [];
ncf_pts = 0;
cf_start=531.25;
cf_end = 538.04;
ifcfplot = 0;
nplots = 0;
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    if (tstamps(1)>tmax)
       break
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
         cp = cp(:);
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
         
         figure(h1)      
%         pvar = plot(x(:)/Chord,dtmp_v(:), 'b.', 'MarkerSize', 10);
%         pvar = plot(xsort/Chord,cp, 'b.', 'MarkerSize', 10);
         pvar = plot(xsort/Chord,cf, 'b.', 'MarkerSize', 10);
%         pvar = plot(xsort,ysort, 'b.', 'MarkerSize', 10);
        
%         set(gca,'Ydir', 'reverse')
%         ylim([-1.1 1.1]);
         xlim ([0. 1]);
%         xlim([-0.01 1.000])
         ylim([-0.004 0.004])    
         grid on   
         hold on
         lgs{1} =  ['T=' num2str(tstamps(it))]; 
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
     



