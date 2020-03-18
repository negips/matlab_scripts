% Get the growth rate from fsi_io

clear
clc
close all

dstn2='/home/prabal/workstation/git_kth/fsi-paper/JFM-latex/imgs';

ifsave = 0;
%svname1='cossu_evol.eps';
%svname2='cossu_peaks.eps';
%svname3='cossu_growth.eps';
%svname4='cossu_omega.eps';

%svname1='lu_evol.eps';
%svname2='lu_peaks.eps';
%svname3='lu_growth.eps';
%svname4='lu_omega.eps';

%svname1='re23_evol.eps';
%svname2='re23_peaks.eps';
%svname3='re23_growth.eps';
%svname4='re23_omega.eps';


outpos = [0.10 0.30 0.40 0.60];

axfs=24;
set(groot,'DefaultAxesFontSize',axfs)
disp(['Default Axes font  size ' num2str(axfs)])

lafs = 30;
lgfs = 24;
mkrsz = 8;
destn = '/home/prabal/workstation/phd_presentations/stability/fsi_linearization/imgs2/';
dstn2 = '/home/prabal/workstation/git_kth/fsi-paper/JFM-latex/imgs/';

fname1 = 'fsi_io.out';
fname2 = '/scratch/negi/git_repos/fsi/stability/baseflow_solve/run_cossu/big_domain/npert_1/fsi_io.out';
fname3 = '/scratch/negi/git_repos/fsi/stability/baseflow_solve/run_cossu/big_domain/nonlinear2/fsi_io.out';

i=1; % 1
file{i}='fsi_lu_pert.out';
legen{i}='Pert';

i=i+1; % 2
file{i}='fsi_lu_base_pert.out';
legen{i}='Base+Pert';

i=i+1; % 3
file{i}='fsi_lu_base_pert_displ.out';
legen{i}='Base+Pert+Displ';

i=i+1; % 4
file{i}='fsi_lu_nl.out';
legen{i}='NL';

% Re=50
i=i+1; % 5
file{i}='fsi_lu_re45_pert.out';
%legen{i}='Re50; Pert';
legen{i}='Linear';

i=i+1; % 6
file{i}='fsi_lu_re45_base_pert.out';
legen{i}='Re50; Base+Pert';

i=i+1; % 7
file{i}='fsi_lu_re45_nl.out';
%legen{i}='Re50; NL';
legen{i}='Nonlinear';

% After pressure bug change
i=i+1; % 8
file{i}='fsi_lu_re45_pert_2.out';
legen{i}='Re50; Pert 2';

i=i+1; % 9
file{i}='fsi_lu_re45_base_pert_2.out';
legen{i}='Re50; Base+Pert 2';

% Non-linear amplitudes

i=i+1; % 10
file{i}='fsi_lu_re45_nl_amp2.out';
legen{i}='Re50; NL; AMP2';

i=i+1; % 11
file{i}='fsi_lu_re45_nl_amp3.out';
legen{i}='Re50; NL; AMP3';

i=i+1; % 12
file{i}='fsi_lu_re45_nl_amp4.out';
legen{i}='Re50; NL; AMP4';

% Including grad(U).eta terms in forces
i=i+1; % 13
file{i}='fsi_lu_re45_all_npert2.out';
legen{i}='Re50; All; npert=2';

i=i+1; % 14
file{i}='fsi_lu_re45_all_npert1.out';
legen{i}='Re50; All; npert=1';

i=i+1; % 15
file{i}='fsi_cossu_re23.out';
legen{i}='Linear'; % Re23.512

i=i+1; % 16
file{i}='fsi_cossu_re23_nl.out';
legen{i}='Nonlinear';% Re23.512

i=i+1; % 17
file{i}='fsi_lu_re45_ugis_pert.out';
legen{i}='Ugis Re45; Lin';

i=i+1; % 18
file{i}='fsi_lu_re45_ugis_NL.out';
legen{i}='Ugis Re45; NL';

i=i+1; % 19
file{i}='fsi_cossu23_n7.out';
legen{i}='Linear'; % Re23.512 n=1/7

i=i+1; % 20
file{i}='fsi_cossu23_n7_nl.out';
legen{i}='Nonlinear';% Re23.512 n=1/7

i=i+1; % 21
file{i}='fsi_ugis_re45_pert.out';
legen{i}='Linear';% Re23.512 n=1/7

i=i+1; % 22
file{i}='fsi_ugis_re45_nl.out';
legen{i}='Nonlinear';% Re23.512 n=1/7

i=i+1; % 23
file{i}='fsi_23_lin.out';
legen{i}='Linear';% Re23.512 n=1/7

i=i+1; % 24
file{i}='fsi_23_nl.out';
legen{i}='Nonlinear';% Re23.512 n=1/7


ind=[21,22];

file=file(ind);
legen=legen(ind);

cols = ['k','k','k','m','c','g','y'];
linst= {'-','none','-','-','-','-','-'};
mkr= {'none','o','.','s','x','*','pentagram'};
mkr_step=1000;


iskip=10;          % no of initial peaks to skip
eskip=0;          % no of end peaks to skip
tstart=-700;
tend  = 800;

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
  ts(i) = plot(time,eta1, 'LineWidth', 2, 'Color', cols(i), 'LineStyle', linst{i}, 'Marker', mkr{i}, 'MarkerSize', mkrsz); hold on
  xlabel('Time', 'FontSize', lafs)
  ylabel('$\eta$', 'FontSize', lafs+4)
  mind = get(ts(i),'MarkerIndices');
  mind2 = mind(1):mkr_step:mind(end);
  set(ts(i),'MarkerIndices',mind2);

  xlim([0 65]);
%  ylim([-2.0 2.0]*10^-6)
  
  [pks locs] = findpeaks(eta);
% Skip negative peaks  
  pks  =  pks(1:2:end);
  locs = locs(1:2:end);
  
  l1 = length(locs);
  locs2 = locs(iskip+1:l1-eskip);

  figure(2)  
  set(gcf,'Units','Normalized')
  set(gcf,'OuterPosition',outpos);
  set(gcf,'Renderer','painters');
  pks2 = pks(iskip+1:l1-eskip);
  pks_time2 = time(locs2);
  ts_pks(i)=semilogy(pks_time2,pks2, 'LineWidth', 2, 'Color', cols(i), 'LineStyle', linst{i}, 'Marker', mkr{i}, 'MarkerSize', mkrsz); hold on
  mind = get(ts_pks(i),'MarkerIndices');
  mind2 = mind(1):4:mind(end);
  set(ts_pks(i),'MarkerIndices',mind2);
 
  xlabel('Time', 'FontSize', lafs)
  ylabel('$\eta^{pks}$', 'FontSize', lafs+4)
 
  
  tosc = diff(pks_time2);
  omg = 2*pi./tosc;
  Omega=mean(omg);
  
  growth = log(pks2(2:end)./pks2(1:end-1))./tosc;
  %growth = diff(pks2)./tosc./pks2(1:end-1);
  Growth=mean(growth);

  growth2 = log(pks2(2:end)./pks2(1:end-1));
 
  time_growth = pks_time2(2:end); 
  figure(3)
  set(gcf,'Units','Normalized')
  set(gcf,'OuterPosition',outpos);
  set(gcf,'Renderer','painters'); 
  plot(time_growth,growth,'LineWidth', 1, 'Color', cols(i), 'LineStyle', linst{i}, 'Marker', mkr{i}, 'MarkerSize', mkrsz); hold on
  ylabel('Growth rate')
  xlabel('Time')

%  xlim([0 1000]);
%  ylim([-14 6]*10^-3)
 

  time_osc = pks_time2(2:end); 
  figure(4)
  set(gcf,'Units','Normalized')
  set(gcf,'OuterPosition',outpos);
  set(gcf,'Renderer','painters');
  plot(time_osc,omg,'LineWidth', 1, 'Color', cols(i), 'LineStyle', linst{i}, 'Marker', mkr{i}, 'MarkerSize', mkrsz); hold on
  ylabel('Angular frequency')
  xlabel('Time')

%  xlim([0 1000]);

  disp(fname)
  disp(['Mean Angular frequency:', num2str(Omega,10)])
  disp(['Mean Growth rate      :', num2str(Growth,10)])

  growth_all{i} = growth;

end

figure(1)
legend(ts,legen,'FontSize',lgfs,'Location','Best')
figure(2)
legend(legen,'FontSize',lgfs,'Location','Best')
figure(3)
legend(legen,'FontSize',lgfs,'Location','Best')
figure(4)
legend(legen,'FontSize',lgfs,'Location','Best')

if (ifsave)
  figure(1)
  SaveFig(gcf,svname1,dstn2,0)

  figure(2)
  SaveFig(gcf,svname2,dstn2,0)

  figure(3)
  SaveFig(gcf,svname3,dstn2,0)

  figure(4)
  SaveFig(gcf,svname4,dstn2,0)
end

