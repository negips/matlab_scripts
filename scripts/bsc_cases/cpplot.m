% Plot cp files

clear
clc
close all


fol = 'beskow_p6';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 0.0;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 1 1])

mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating


nplots = 0;
for i = 1:1
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);
    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast)
         if nplots>0
           delete(pvar)
         end   
         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);
         dtmp_v = sdata(3).data(:,:,it);    
         pvar = plot(x(:),dtmp_v(:), 'b.', 'MarkerSize', 6);
         set(gca,'Ydir', 'reverse')
         ylim([-3.5 1.1])
         hold on
         lgs{1} =  ['T=' num2str(tstamps(it))]; 
         lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'NorthEast', 'Fontsize', lfs);
         if nplots == 0   
           ylabel('$C_{p}$', 'Interpreter', 'Latex', 'Fontsize', fs);
           xlabel('$x/C$', 'Interpreter', 'Latex', 'Fontsize', fs);
         end   
         nplots = nplots+1;   
         mov(nplots) = getframe(gcf);
 
      end
      tlast = tstamps(it);
      pause(0.001)
    end
  end
end           

mov2 = mov(1:nplots);      
movie2avi(mov2, 'cp_movie.avi', 'compression', 'None', 'fps', 20, 'quality', 50);
%      for ies=1:selt
%      %     scatter(x(:,i),sdata(1).data(:,i,1), '.')
%      %     scatter(x(:,i),sdata(1).data(:,i,maxtsaves), 'd')
%           scatter(sdata(1).data(:,ies,maxtsaves),sdata(3).data(:,ies,maxtsaves), 'd')
%      
%      end
%      set(gca,'Ydir','reverse')
     



