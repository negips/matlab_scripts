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
tlast = 19.1;
tmax = 20.50001;
maxframes = nfiles*100;

% Load surface normals
snormals = importdata('surf_normals.67');

x_imp = snormals.data(:,1);
[x_imp I] = sort(x_imp);
y_imp = snormals.data(I,2);

snx = -snormals.data(I,4);
sny = -snormals.data(I,5);

topind = sny>0;
snx_top = snx(find(topind));
sny_top = sny(find(topind));
xt_imp  = x_imp(find(topind));
yt_imp  = y_imp(find(topind));

snx_bot = snx(find(~topind));
sny_bot = sny(find(~topind));
xb_imp  = x_imp(find(~topind));
yb_imp  = y_imp(find(~topind));

for i=1:length(snx_top)
  stx_top(i) = sny_top(i);
  sty_top(i) = -snx_top(i);
end

for i=1:length(snx_bot)
  stx_bot(i) = -sny_bot(i);
  sty_bot(i) = snx_bot(i);
end



h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

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
tall  = [];

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

         if nplots>0
           delete(pvar)
         end
         dtmpx = sdata(1).data(:,:,it);
         dtmpy = sdata(2).data(:,:,it);
         dtmp_v = sdata(3).data(:,:,it);    
%         dtmp_v = -sdata(5).data(:,:,it);    
         cf = -sdata(5).data(:,:,it);

         xnew=dtmpx(:);
         ynew=dtmpy(:);
         cf = cf(:);    

%        Reallign and subsample
         if (mod(nplots,5)==0)   
           arg = omg*(tstamps(it)-pitch_start);
           theta = -pitch_amp*sin(arg);    % negative angle to go backwards
%          Clockwise rotation by theta 
           rot = [cos(theta) sin(theta);
                  -sin(theta) cos(theta)];   
           axis_x0=0.35;
           axis_y0=0.034;            
           coords = rot*[dtmpx(:)'-axis_x0; dtmpy(:)'-axis_y0];
           xnew = coords(1,:)+axis_x0;
           ynew = coords(2,:)+axis_y0; 
 
           [xsort ind] = sort(xnew);
           xsort = xsort'; 
           ysort = ynew(ind)';
           cfnew = cf;
           cf = cfnew(ind);


%          Truncate space to top surface
           ind2=find(topind);   
 
           cfall = [cfall cf(ind2)];
           tall = [tall tstamps(it)];
 
           xnew = xsort;
           ynew = ysort;
         end    
%%
         figure(h1)      
%         pvar = plot(x(:)/Chord,dtmp_v(:), 'b.', 'MarkerSize', 10);
%         pvar = plot(xsort/Chord,cf, 'b.', 'MarkerSize', 15);
         pvar = plot(xnew,cf, 'b.', 'MarkerSize', 15);
%         set(gca,'Ydir', 'reverse')
%         ylim([-3.5 1.1]);
%         xlim([-0.01 1.0])
%         ylim([-0.003 0.0040])    
         grid on   
         hold on
         lgs{1} =  ['T=' num2str(tstamps(it))]; 
         lg = legend(pvar,lgs, 'FontSize', lfs, 'Location', 'North', 'Fontsize', lfs, 'Box', 'off');


         if nplots == 0 
           ylabel('\tau_{w}', 'Interpreter', 'tex', 'Fontsize', fs);
           xlabel('x/C', 'Interpreter', 'tex', 'Fontsize', fs);
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


[lg nkryl]=size(cfall);

%nkryl = 30;

[U,S,W]=svd(cfall(:,1:nkryl-1),0);
stilde=(U')*(cfall(:,2:nkryl))*W*inv(S);
[evec ritz] = eig(stilde);
ev = diag(ritz);
ritz_r=real(ev);
ritz_i=imag(ev);

evec_r = real(evec);
evec_i = imag(evec);

deltaT=mean(diff(tall));      % should be constant

eigen_r = log(sqrt(ritz_r.^2 + ritz_i.^2))/deltaT;
eigen_i = atan2(ritz_i,ritz_r)/deltaT;
figure(10)
plot(eigen_r,eigen_i, '.b', 'MarkerSize', 16); hold on
grid on

dmd_r=zeros(lg,nkryl-1);
dmd_i=zeros(lg,nkryl-1);
for i=1:nkryl-1
  dmd_r(:,i)=U*evec_r(:,i);
  dmd_i(:,i)=U*evec_i(:,i);
end

[vals inds]=sort(eigen_r,1,'descend');

for j=1:length(inds)
  i=inds(j);
  figure(11); 
  plot(xsort(find(topind)),dmd_r(:,i), '.k', 'MarkerSize', 16);
  legend(num2str(j), 'Location', 'SouthWest')
  [eigen_r(i) eigen_i(i)]
  pause(3)
end




