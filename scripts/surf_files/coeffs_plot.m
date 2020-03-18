% Plot cp files

clear
clc
close all

% addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/'
% addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

fol = 'saab_sfd_k1.0';
ifhdr = 1;
fs = 16;                % fontsize
lfs = 16;               % legend fontsize
lafs = 24;              % Latex font size
ifcols = 1;
destn = 'plots/';

[sfiles tout] = LoadSurfFiles(fol);

nfiles = length(sfiles);
tlast = 60.00;
tmax = 500.00001;
maxframes = nfiles*100;

h1=figure('units','normalized','outerposition',[0 0 0.4 0.6]);

dragxp = [];
dragyp = [];

dragxv = [];
dragyv = [];

time   = [];

cfavg = [];
cfavgx = [];
cfavgy = [];
ncf_pts = 0;
cf_start=531.25;
cf_end = 538.04;
ifcfplot = 1;
nplots = 0;
for i = 1:nfiles
  if (tout(i)>=tlast)
    fname = sfiles{i};
       
    ndim=2; 
    [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr,ndim);

    dragxp = [dragxp; sintegrals(:,1,3)]; 
    dragyp = [dragyp; sintegrals(:,2,3)]; 

    dragxv = [dragxv; sintegrals(:,1,5)]; 
    dragyv = [dragyv; sintegrals(:,2,6)];

    time   = [time; tstamps]; 

    K=1.0;
    omg=K*2;
    alpha   = 0.8*sin(omg*time);

    phi=-2.40;
    alpha_e = 0.8*sin(omg*time+phi);

    plot(time,dragyp);hold on
    plot(time,alpha, '--k');hold on
    plot(time,alpha_e, '--r');hold off
    grid on
    pause(0.01)
  end
end           

     



