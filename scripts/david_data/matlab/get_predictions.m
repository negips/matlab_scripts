% clear
% clc
% close all

addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

destn = 'plots/';
ifcols= 1;

fs=16;
lfs=12;

all = load('all_predictions.mat');;

[k ind] = sort(all.kall2);
k = [0 k];
phi = [0 all.phiall2(ind)];
gamma = [0 all.gammaall2(ind)];
intgbydalpha = [0 all.intgbydalpha2(ind)];
intgnorm = [0 all.intgnorm2(ind)];

nsmooths=50; 
movk=k;
movphi=phi;
movgamma=gamma;
movintg=intgbydalpha;
movintgnorm=intgnorm;
for j=1:nsmooths
  order=1;
  framelen=3;  
  movk = sgolayfilt(movk,order,framelen);
  movphi = sgolayfilt(movphi,order,framelen);
  movgamma = sgolayfilt(movgamma,order,framelen);
  movintg = sgolayfilt(movintg,order,framelen);
  movintgnorm = sgolayfilt(movintgnorm,order,framelen);
end

figure(10)
plot(k,phi*180/pi, '*'); hold on
plot(movk,movphi*180/pi, '-k', 'LineWidth', 2)
plot(movk,1.2*movphi*180/pi, '--r', 'LineWidth', 1)
plot(movk,0.8*movphi*180/pi, '--r', 'LineWidth', 1)
ylabel('\phi', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
filename=['model_k-phi.eps'];
SaveFig(gcf,filename, destn, ifcols)

figure(11)
plot(k,gamma*180/pi, '*'); hold on
plot(movk,movgamma*180/pi, '-k', 'LineWidth', 2)
plot(movk,1.3*movgamma*180/pi, '--r', 'LineWidth', 1)
plot(movk,0.7*movgamma*180/pi, '--r', 'LineWidth', 1)
ylabel('\phi', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)

figure(12)
plot(k,intgbydalpha, 'd'); hold on
plot(movk,movintg, '-k', 'LineWidth', 2)
plot(movk,1.15*movintg, '--r', 'LineWidth', 1)
plot(movk,0.85*movintg, '--r', 'LineWidth', 1)
ylabel('Integ const', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
filename=['model_k-intg.eps'];
SaveFig(gcf,filename, destn, ifcols)

figure(13)
plot(k,intgnorm, 'd'); hold on
plot(movk,movintgnorm, '-k', 'LineWidth', 2)
plot(movk,1.1*movintgnorm, '--r', 'LineWidth', 1)
plot(movk,0.9*movintgnorm, '--r', 'LineWidth', 1)
ylabel('Normalized Integ const', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
filename=['model_k-intgnorm.eps'];
SaveFig(gcf,filename, destn, ifcols)

