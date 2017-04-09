% Plot cp files

clear
clc
close all

% addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 're750k_aoa44';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 4.50;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

nplots = 0;
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);
    xmin=min(x(:));
    xmax=max(x(:));
    Chord=xmax-xmin;
    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast)
         if nplots>0
           delete(pvar)
         end   
         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);
%         dtmp_v = sdata(3).data(:,:,it);    
         dtmp_v = -sdata(5).data(:,:,it);    

         pvar = plot(x(:)/Chord,dtmp_v(:), 'b.', 'MarkerSize', 6);
%         set(gca,'Ydir', 'reverse')
%         ylim([-3.5 1.1]);
         hold on
         lgs{1} =  ['T=' num2str(tstamps(it))]; 
         lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'North', 'Fontsize', lfs, 'Box', 'off');
         if nplots == 0 
           ylabel('C_{p}', 'Interpreter', 'tex', 'Fontsize', fs);
           xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
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

%mov2 = mov(1:nplots); 
%movie2avi(mov2, 'cp_movie.avi', 'compression', 'None', 'fps', 15, 'quality', 75);
%      for ies=1:selt
%      %     scatter(x(:,i),sdata(1).data(:,i,1), '.')
%      %     scatter(x(:,i),sdata(1).data(:,i,maxtsaves), 'd')
%           scatter(sdata(1).data(:,ies,maxtsaves),sdata(3).data(:,ies,maxtsaves), 'd')
%      
%      end
%      set(gca,'Ydir','reverse')
     



