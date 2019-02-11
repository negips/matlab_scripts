% Get the growth rate from fsi_io

clear
clc
close all

lafs = 20;
lgfs = 16;

destn = '/home/prabal/workstation/phd_presentations/stability/fsi_linearization/imgs2/';

fname1 = 'fsi_io.out';
fname2 = '/scratch/negi/git_repos/fsi/stability/baseflow_solve/run_cossu/big_domain/npert_1/fsi_io.out';
fname3 = '/scratch/negi/git_repos/fsi/stability/baseflow_solve/run_cossu/big_domain/nonlinear2/fsi_io.out';

i=1; % 1
file{i}='eigenvalues_cossuN7.txt';
legen{i}='Cossu; N=7';

i=i+1; % 2
file{i}='eigenvalues_cossuN9.txt';
legen{i}='Cossu; N=9';

i=i+1; % 3
file{i}='eigenvalues_lurotatedN7.txt';
legen{i}='Rotated Splitter-plate; N=7';

i=i+1; % 4
file{i}='eigenvalues_lurotatedN9.txt';
legen{i}='Rotated Splitter-plate; N=9';

i=i+1; % 5
file{i}='eigenvalues_ugis45.txt';
legen{i}='Ugis Re45; N=7';

i=i+1; % 6
file{i}='eigenvalues_ugis156.txt';
legen{i}='Ugis Re156; N=11';

i=i+1; % 7
file{i}='eigenvalues_lu100N11.txt';
legen{i}='Lu (2016) Re100; N=11';

ind=[7];

file=file(ind);
legen=legen(ind);

cols = ['b','r','k','m','c','g','y'];
linst= {'-','--','-','-','-','-','-'};
mkrs = {'o','*','s','d','.'};

nfiles=length(ind);

figure(1)
set(gcf,'Units','Normalized')
figpos = [0.35 035 0.45 0.60];
set(gcf,'OuterPosition', figpos);
set(gcf,'Renderer','Painters')
for i=1:nfiles

  fname = file{i};
  egv = importdata(fname);
  ritzr   = egv.data(:,3);
  ritzi   = egv.data(:,2);
  eigenr  = egv.data(:,5);
  eigeni  = egv.data(:,4);

  ind = find(eigenr>=-1e-12);

  spec(i) = plot(eigenr(ind),eigeni(ind),'LineWidth', 2, 'Color', cols(i), 'LineStyle', 'none', 'Marker', mkrs{i}); hold on
  xlabel('$\omega_{r}$', 'FontSize', lafs)
  ylabel('$\omega_{i}$', 'FontSize', lafs)

end

figure(1)

xlims = get(gca,'XLim');
zeroline = plot(xlims,[0 0], '--', 'LineWidth', 2, 'Color', 'k');

legend(spec,legen,'FontSize',lgfs,'Location','Best')
xlim(xlims)

SaveFig(gcf,'lu100_spectrum',destn,1)
