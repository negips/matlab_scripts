% Create a model for the static Cm/Cz curves

clear
clc
close all

%addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

base = '../';
%base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
% fname = 'u24_a3_d14_fsteps.h5';
lfol = length(fol);
destn = 'plots/';
ifcols= 1;

fs = 16;    % fontsize
lfs = 12;   % legend font size
lw = 1;     % linewidth

c=0.5;
nu = 1.568E-5;

%fname='u30_d14_alphasweep1.h5';
fname='u24_d14_alphasweep.h5';
%fname='u30_d8_alphasweep.h5';
%fname='u24_d8_alphasweep.h5';
uoo = 24;
deltacase=14;
Re=uoo*c/nu;
if (Re<600000)
  re_leg = '565k';
elseif (Re>600000 && Re < 900000)
  re_leg = '765k';
else
  re_leg = '950k';
end        

hfile = [folder fname];

[segments] = split_segments(hfile, uoo, deltacase);
allcount = 1;
legs{allcount} = ['Re=' re_leg '; \delta=', num2str(deltacase)];

cz_all = [];
cm_all = [];
alpha_all = [];
time_all = [];

alpha_min = 0;         % degrees
alpha_max = 10;         % degrees

nsegs = length(segments);
for iseg=1:nsegs
  qtime=segments(iseg).qtime;
  alpha=segments(iseg).alpha*180/pi;
  p_time = segments(iseg).ptime;
  p_cz = segments(iseg).Cz;
  p_cm = segments(iseg).Cm;

%  figure(1)
%  plot_aoa(allcount) = plot(qtime,alpha, '.', 'Color', 'b'); 
%  hold on

  ind1=alpha>alpha_min;
  ind2=alpha<alpha_max;
  ind3=find(ind1.*ind2);
  alpha2=alpha(ind3);
  qtime2=qtime(ind3);

  q_cz = interp1(p_time,p_cz,qtime2,'pchip');
  q_cm = interp1(p_time,p_cm,qtime2,'pchip');

  cz_all = [cz_all q_cz'];
  cm_all = [cm_all q_cm'];
  alpha_all = [alpha_all alpha2'];
  time_all = [time_all qtime2'];

end

[alpha3 ind4] = sort(alpha_all);
cm3 = cm_all(ind4);
cz3 = cz_all(ind4);

nsmooths=50; 
movalpha=alpha3;
movcm=cm3;
movcz=cz3;
for j=1:nsmooths
  order=1;
  framelen=201;  
  movalpha = sgolayfilt(movalpha,order,framelen);
  movcm = sgolayfilt(movcm,order,framelen);
  movcz = sgolayfilt(movcz,order,framelen);
end

figure(3)
plot_cq1 = plot(alpha3,cm3, '.', 'Color', 'b'); hold on;
%plot_cq = plot(movalpha,movcm, '--', 'Color', 'm', 'LineWidth',3);
ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
hold on
xlim([-1 11])

figure(4)
plot_cq1 = plot(alpha3,cz3, '.', 'Color', 'b'); hold on;
%plot_cq = plot(movalpha,movcz, '--', 'Color', 'm', 'LineWidth',3);
ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
hold on
xlim([-1 11])

%% Reduce resolution and save
npts=10000;
alpha = linspace(alpha_min,alpha_max,npts);
cm = interp1(movalpha,movcm,alpha,'pship');
cz = interp1(movalpha,movcz,alpha,'pship');

figure(3)
plot_cq2=plot(alpha,cm,'-d', 'Color', 'c', 'LineWidth', 3);
filename=['static_model_cm.eps'];
filename = [re_leg '_' filename];
SaveFig(gcf,filename, destn, ifcols)

figure(4)
plot_cq2=plot(alpha,cz,'-d', 'Color', 'c', 'LineWidth', 3);
filename=['static_model_cz.eps'];
filename = [re_leg '_' filename];
SaveFig(gcf,filename, destn, ifcols)


% save([num2str(deltacase) '_static_models_' num2str(re_leg) '.mat'], 'alpha', 'cm', 'cz', 'deltacase', 'uoo', 'base', 'fol', 'folder', 'fname')






