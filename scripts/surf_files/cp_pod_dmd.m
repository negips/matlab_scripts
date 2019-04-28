% Plot cp files

clear
clc
%close all

mkr='.';
mc ='k';
msiz=16;
nskip=4;

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 'naca0012_re77k_large';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
ifcols = 1;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 50.0;
tmax = 200.50001;
maxframes = nfiles*100;

%h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

%mov(1:maxframes) = struct('cdata', [],'colormap', []);            % Just allocating
%mov = VideoWriter('cp_movie.avi');

nplots = 0;

kred=0.5;
Uinf=1.0;
semichord=0.5;
omg = kred*Uinf/semichord;
pitch_start=0;
pitch_amp=1.3*pi/180;

cfall = [];
prall = [];
tall  = [];

for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
        
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr);

    if (tstamps(1)>tmax)
       break
    end         

    xmin=min(x(:));
    xmax=max(x(:));
    Chord=xmax-xmin;
%    Chord=1;

%   Very first x/y - use as ref  
    if nplots==0
      xref = x;
      yref = y;
    end    

    for it = 1:length(tstamps)
      if (tstamps(it)>=tlast)
         
         if (tstamps(it)>tmax)
            break
         end

         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);
         dtmp_v = sdata(3).data(:,:,it);    
%         dtmp_v = -sdata(5).data(:,:,it);    
         cf = -sdata(5).data(:,:,it);
         pr = sdata(3).data(:,:,it); 

         xnew=dtmpx(:);
         ynew=dtmpy(:);
         cf = cf(:);
         pr = pr(:);

%        Reallign and subsample
         if (mod(nplots,nskip)==0)   
 
%           xsort = xsort'; 
%           ysort = ynew(ind)';
%           cfnew = cf;
%           cf = cfnew(ind);
 
           cfall = [cfall cf];
           prall = [prall pr];
           tall = [tall tstamps(it)];

         end    
%%

         svfname = sprintf('%0.4d', nplots);   
         svfname = ['cf_fig' svfname '.eps'];
         destn = 'plots/';   
%         SaveFig(gcf, svfname, destn, 1)
      end
      
      nplots = nplots+1;
      tlast = tstamps(it);
%      pause(0.01)
    end           % it
  end
end           


V = prall;

[lg nkryl]=size(V);

%nkryl = 30;

[U,S,W]=svd(V(:,1:nkryl-1),0);
s=diag(S);
stilde=(U')*(V(:,2:nkryl))*W*inv(S);
[evec ritz] = eig(stilde);
ev = diag(ritz);
ritz_r=real(ev);
ritz_i=imag(ev);

figure(8)
semilogy(s, mc, 'Marker', '.', 'MarkerSize', 10); hold on
title('Singular Values')

figure(9)
plot(ritz_r,ritz_i, [mkr mc], 'MarkerSize', msiz); hold on
grid on


evec_r = real(evec);
evec_i = imag(evec);

deltaT=mean(diff(tall));      % should be constant

eigen_r = log10(sqrt(ritz_r.^2 + ritz_i.^2))/deltaT;
eigen_i = atan2(ritz_i,ritz_r)/deltaT;
figure(10)
plot(eigen_r,eigen_i, [mkr mc], 'MarkerSize', msiz); hold on
grid on

dmd_r=zeros(lg,nkryl-1);
dmd_i=zeros(lg,nkryl-1);
for i=1:nkryl-1
  dmd_r(:,i)=U*evec_r(:,i);
  dmd_i(:,i)=U*evec_i(:,i);
end

%[vals inds]=sort(eigen_r,1,'descend');

%for j=1:length(inds)
%  i=inds(j);
%  figure(11); 
%%  plot(xsort(find(topind)),dmd_r(:,i), '.k', 'MarkerSize', 16);
%  plot(xnew,dmd_r(:,i), '.b', 'MarkerSize', 16); hold on
%  plot(xnew,dmd_i(:,i), '.r', 'MarkerSize', 16); hold off
%
%
%  legend(num2str(j), 'Location', 'SouthWest')
%  [i eigen_r(i) eigen_i(i)]
%  pause(2)
%end

figure(12)
plot(xnew,V(:,1), '.k', 'MarkerSize',16); hold on
i=1;
scale=mean(W(:,i))*s(i);
plot(xnew,U(:,i)*scale, '.b', 'MarkerSize',16)

figure(13)
plot(tall(1:end-1),W(:,1), 'b'); hold on
%plot(tall(1:end-1),W(:,2), 'r'); hold on








