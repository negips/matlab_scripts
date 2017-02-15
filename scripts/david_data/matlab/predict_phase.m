% Trying to predict phase plots from previous data

clear
clc
close all

all = load('all_predictions.mat');;

[k ind] = sort(all.kall2);
k = [0 k];
phi = [0 all.phiall2(ind)];
intgbydalpha = [0 all.intgbydalpha2(ind)];

nsmooths=50; 
movk=k;
movphi=phi;
movintg=intgbydalpha;
for j=1:nsmooths
  order=1;
  framelen=3;  
  movk = sgolayfilt(movk,order,framelen);
  movphi = sgolayfilt(movphi,order,framelen);
  movintg = sgolayfilt(movintg,order,framelen);
end

figure(10)
plot(k,phi*180/pi, '*'); hold on
plot(movk,movphi*180/pi, 'LineWidth', 2)

figure(11)
plot(k,intgbydalpha, 'd'); hold on
plot(movk,movintg, 'LineWidth', 2)

static_model = load('14_static_models_765k.mat');

mean_aoa=2.9;
dalpha=0.9;
k=0.3;           % min value is 0.02

omegat = linspace(0,4*pi,1000);
intg_const = interp1(movk,movintg,k,'linear');
intg_const = 1.0*intg_const*dalpha*pi/180;
philag = interp1(movk,movphi,k,'linear');

[intg_const philag*180/pi]

aoa = mean_aoa + dalpha*sin(omegat);
pitch = intg_const*cos(omegat);

alphalagg =  mean_aoa + dalpha*sin(omegat + philag);
cz_lag = interp1(static_model.alpha,static_model.cz,alphalagg,'linear');

cz = pitch + cz_lag;

figure(12)
plot(aoa,cz, 'LineWidth', 2); hold on
plot(static_model.alpha,static_model.cz, '--k', 'LineWidth', 2)
xmin = mean_aoa-dalpha*1.2;
xmax = mean_aoa+dalpha*1.2;
xlim([xmin xmax])
title('Phase Potrait')
grid on

figure(13)
plot(omegat,cz, 'LineWidth', 2);
ax1=gca;
%ax1_pos=get(ax1,'Position');
%ax2=axes('Position',ax1_pos, 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'Color', 'none')
%plot(omegat,aoa, 'Parent',ax2,'Color','r', 'LineStyle','--','LineWidth',2);
%xmin = mean_aoa-dalpha*1.2;
%xmax = mean_aoa+dalpha*1.2;
%xlim([xmin xmax])
title('Cz-time')






