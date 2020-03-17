% Plot cp files

clear
clc
close all

% addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 'seck0.4_t0.25/';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
ifplot = 0;             % plot individual wall profiles
tlast = 1.9635;          % start from this time
tstart0 = tlast;
tend = 3.2385;             % stop at this time
destn = 'plots/';
ifcp = 1;               % plot pressure instead of cf
ifdatasave=1;           % save data into a mat file
  datafile='sec_k04_t025.mat';
lafs = 22;              % Latex font size

U0=1.;
kred=0.0;
chord=1.0;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=1.; %2*pi/omega;
%Tosc=1;
ptch_amp = 0.0;
ptch_start = 1.9635;
axis_x0 = 0.35;
axis_y0 = 0.034;
phase=0.;% -pi/2;

pitch_end=4.5*Tosc + ptch_start;

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
maxframes = nfiles*100;

if (ifplot)
  h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);
end

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

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
  if (tout(i)>=tlast && tout(i)<=tend) 
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast && tstamps(it)<=tend)
         if ifplot && nplots>0
           delete(pvar)
         end

         xx = sdata(1).data(:);
         yy = sdata(2).data(:);
         [x1,mini] = min(xx(:));
         [x2,maxi] = max(xx(:));
         Chord=1.0;
         y1 = yy(mini);
         y2 = yy(maxi);

         slope = (y2-y1)/(x2-x1);

         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);
%         dtmp_v = sdata(3).data(:,:,it);
         dtmp_cp = sdata(3).data(:,:,it);
         pmax = max(dtmp_cp(:));

         dtmp_v = -sdata(5).data(:,:,it);
         dtmp_cfx = -sdata(5).data(:,:,it);
         dtmp_cfy = -sdata(6).data(:,:,it);
      
         [xsort ind1] = sort(dtmpx(:));
         ysort = dtmpy(ind1);

         yline = y1 + slope*(xsort-x1);

         ind = find(ysort>=yline);  % These points are above the line

         xsort = xsort(ind);
         ysort = ysort(ind);

         cfx = dtmp_cfx(ind1);
         cfy = dtmp_cfy(ind1); 
         cp = dtmp_cp(ind1);

         cfx = cfx(ind);
         cfy = cfy(ind);
         cp  = cp(ind);

         cf = cfx;
                         
         x_min = x1;
         x_max = x2;

         if (ifplot)
           if (ifcp)
             pvar = plot(xsort,cp, 'b.', 'MarkerSize', 6);
             grid on
             set(gca,'Ydir', 'reverse')
%             ylim([-1.1 1.1]);
             ylabel('C_{p}', 'Interpreter', 'tex', 'Fontsize', fs);
             xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
           else        
             pvar = plot(xsort,cfx./max(cfx), 'b.', 'MarkerSize', 8);
%             pvar = semilogy(xsort,(abs(cfx)+1.0e-12), 'b.', 'MarkerSize', 8);
             grid on
%             ylim([-0.0035 0.005]);
             ylim([-1.5 1.5])
             xlim([x_min x_max])
%             if nplots == 0 
               ylabel('C_{f}', 'Interpreter', 'tex', 'Fontsize', fs);
               xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
%             end   
           end 
           lgs{1} =  ['T=' num2str(tstamps(it))]; 
           lg = legend(pvar,lgs, 'FontSize', lafs, 'Location', 'SouthWest', 'Box', 'off');

         end      % ifplot

         tstamps(it)=tstamps(it);

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
  save(datafile)
  disp(['Data Saved to ' datafile])
end

%splot_impulse


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
     



