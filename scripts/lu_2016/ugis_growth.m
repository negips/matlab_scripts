% Get the growth rate from fsi_io

clear
clc
close all

dstn2='/home/prabal/workstation/git_kth/fsi-paper/JFM-latex/imgs';

ifsave = 1;
%svname1='cossu_evol.eps';
%svname2='cossu_peaks.eps';
%svname3='cossu_growth.eps';
%svname4='cossu_omega.eps';

svname1='ugis_log.eps';

outpos = [0.10 0.30 0.40 0.60];

axfs=24;
set(groot,'DefaultAxesFontSize',axfs)
disp(['Default Axes font  size ' num2str(axfs)])

lafs = 30;
lgfs = 24;
mkrsz = 8;
destn = '/home/prabal/workstation/phd_presentations/stability/fsi_linearization/imgs2/';
dstn2 = '/home/prabal/workstation/git_kth/fsi-paper/JFM-latex/imgs/';

i=1; % 1
file{i}='fsi_ugis_re45_pert.out';
legen{i}='Linear';% Re45

i=i+1; % 2
file{i}='fsi_ugis_re45_nl.out';
legen{i}='Nonlinear';% Re45

ind=[1,2];

file=file(ind);
legen=legen(ind);

cols = ['k','k','k','m','c','g','y'];
linst= {'-','none','-','-','-','-','-'};
mkr= {'none','o','.','s','x','*','pentagram'};
mkr_step=10000;


iskip=0;          % no of initial peaks to skip
eskip=0;          % no of end peaks to skip
tstart=50;
tend  =350;

%re45
%tstart=950;
%tend=1550;

nfiles=length(ind);
growth_all = [];

for i=1:nfiles

  if ind(i)==4
    eskip=8;
  elseif ind(i)==16
    eskip=53;
  end    

  fname = file{i};
  fsi = importdata(fname);
  time=fsi.data(:,2);
  time=time-min(time);
  eta1 = fsi.data(:,4);
  eta  = abs(fsi.data(:,4)); % + 1e-15;
  etav = fsi.data(:,5);

  if (tstart>0)
    ind1=time>=tstart;
  else
    ind1=ones(length(time),1);
  end  

  if (tend>0)
    ind2=time<=tend;
  else
    ind2=ones(length(time),1);
  end
  ind3=find(ind1.*ind2);
  time=time(ind3);
  eta1=eta1(ind3);
  eta = eta(ind3);
  etav= etav(ind3);
 
  figure(1)
  set(gcf,'Units','Normalized')
  set(gcf,'OuterPosition',outpos);
  set(gcf,'Renderer','painters');
  ts(i) = semilogy(time,abs(eta1), 'LineWidth', 2, 'Color', cols(i), 'LineStyle', linst{i}, 'Marker', mkr{i}, 'MarkerSize', mkrsz); hold on
  xlabel('Time', 'FontSize', lafs)
  ylabel('$\eta$ (radians)', 'FontSize', lafs+4)
  mind = get(ts(i),'MarkerIndices');
  mind2 = mind(1):mkr_step:mind(end);
  set(ts(i),'MarkerIndices',mind2);

%  xlim([0 65]);
%  ylim([-2.0 2.0]*10^-6)
  
%  [pks locs] = findpeaks(eta);
% Skip negative peaks  
%  pks  =  pks(1:2:end);
%  locs = locs(1:2:end);
  
%  l1 = length(locs);
%  locs2 = locs(iskip+1:l1-eskip);
%
%  figure(2)  
%  set(gcf,'Units','Normalized')
%  set(gcf,'OuterPosition',outpos);
%  set(gcf,'Renderer','painters');
%  pks2 = pks(iskip+1:l1-eskip);
%  pks_time2 = time(locs2);
%  ts_pks(i)=semilogy(pks_time2,pks2, 'LineWidth', 2, 'Color', cols(i), 'LineStyle', linst{i}, 'Marker', mkr{i}, 'MarkerSize', mkrsz); hold on
%  mind = get(ts_pks(i),'MarkerIndices');
%  mind2 = mind(1):4:mind(end);
%  set(ts_pks(i),'MarkerIndices',mind2);
% 
%  xlabel('Time', 'FontSize', lafs)
%  ylabel('$\eta^{pks}$', 'FontSize', lafs+4)

end

figure(1)
legend(ts,legen,'FontSize',lgfs,'Location','Best')
%figure(2)
%legend(legen,'FontSize',lgfs,'Location','Best')

if (ifsave)
  figure(1)
  SaveFig(gcf,svname1,dstn2,0)
end

