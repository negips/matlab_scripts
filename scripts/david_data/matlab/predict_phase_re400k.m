% Trying to predict phase plots from previous data

clear
clc
close all

fs=16;
lfs=12;

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
plot(movk,movphi*180/pi, '-k', 'LineWidth', 2)
plot(movk,1.2*movphi*180/pi, '--r', 'LineWidth', 1)
plot(movk,0.8*movphi*180/pi, '--r', 'LineWidth', 1)
ylabel('\phi', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)


figure(11)
plot(k,intgbydalpha, 'd'); hold on
plot(movk,movintg, '-k', 'LineWidth', 2)
plot(movk,1.15*movintg, '--r', 'LineWidth', 1)
plot(movk,0.85*movintg, '--r', 'LineWidth', 1)
ylabel('Integ const', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)

xfoil = importdata('polar_re4e5_ed36f128+14.dat');
%static_model = load('14_static_models_765k.mat');
static_model.alpha=xfoil.data(:,1);
static_model.cz=xfoil.data(:,2);

mean_aoa=4.2;
dalpha=1.5;
k=0.5;           % min value is 0.02
U0=1.;
c=1.0;
semichord=c/2;
omega=k*U0/semichord;

omegat = linspace(0,4*pi,1000);
time=omegat/omega;
intg_const = interp1(movk,movintg,k,'linear', 'extrap');
intg_const = 1.2*intg_const*dalpha*pi/180;
intg_min = intg_const*0.85;
intg_max = intg_const*1.15;

philag = interp1(movk,movphi,k,'linear', 'extrap');
philag_min=0.8*philag;
philag_max=1.2*philag;

[intg_const philag*180/pi]

aoa = mean_aoa + dalpha*sin(omegat);
pitch = intg_const*cos(omegat);
pitch_min = intg_min*cos(omegat);
pitch_max = intg_max*cos(omegat);

alphalagg =  mean_aoa + dalpha*sin(omegat + philag);
alphalagg_min =  mean_aoa + dalpha*sin(omegat + philag_min);
alphalagg_max =  mean_aoa + dalpha*sin(omegat + philag_max);

cz_lag = interp1(static_model.alpha,static_model.cz,alphalagg,'linear', 'extrap');
cz_lagmin = interp1(static_model.alpha,static_model.cz,alphalagg_min,'linear', 'extrap');
cz_lagmax = interp1(static_model.alpha,static_model.cz,alphalagg_max,'linear', 'extrap');

cz = pitch + cz_lag;
cz_1 = pitch_min+cz_lagmin;
cz_2 = pitch_max+cz_lagmax;
cz_3 = pitch_min+cz_lagmax;
cz_4 = pitch_max+cz_lagmin;

col1 = lines(5);
legs_12 = {'Cz' 'Cz1' 'Cz2' 'Cz3' 'Cz4' 'Static Cz'};

figure(12)
plot(aoa,cz, 'LineWidth', 4, 'Color', 'k'); hold on
plot(aoa,cz_1, '-.', 'LineWidth', 2, 'Color', col1(2,:));
plot(aoa,cz_2, '-.', 'LineWidth', 2, 'Color', col1(3,:));
plot(aoa,cz_3, '-.', 'LineWidth', 2, 'Color', col1(4,:));
plot(aoa,cz_4, '-.', 'LineWidth', 2, 'Color', col1(5,:));
plot(static_model.alpha,static_model.cz, '-b', 'LineWidth', 2)
legend(legs_12, 'Interpreter', 'tex', 'FontSize', lfs)
xmin = mean_aoa-dalpha*1.4;
xmax = mean_aoa+dalpha*1.4;
xlim([xmin xmax])
title('Phase Potrait')
ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
grid on

legs_13 = {'Cz' 'Cz1' 'Cz2' 'Cz3' 'Cz4'};
figure(13)
plot(time,cz, 'LineWidth', 3, 'Color', 'k'); hold on
plot(time,cz_1, '-.', 'LineWidth', 2, 'Color', col1(2,:));
plot(time,cz_2, '-.', 'LineWidth', 2, 'Color', col1(3,:));
plot(time,cz_3, '-.', 'LineWidth', 2, 'Color', col1(4,:));
plot(time,cz_4, '-.', 'LineWidth', 2, 'Color', col1(5,:));
legend(legs_13, 'Interpreter', 'tex', 'FontSize', lfs)

ax1=gca;
%ax1_pos=get(ax1,'Position');
%ax2=axes('Position',ax1_pos, 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'Color', 'none')
%plot(omegat,aoa, 'Parent',ax2,'Color','r', 'LineStyle','--','LineWidth',2);
%xmin = mean_aoa-dalpha*1.2;
%xmax = mean_aoa+dalpha*1.2;
%xlim([xmin xmax])
title('Cz-time')
ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('t', 'Interpreter', 'tex', 'FontSize', fs)


%% Xfoil
ifxfoil = 0;
if ifxfoil
  xfoil = importdata('polar_re765k_ed36f128+14.dat');
  figure(12)
  plot(xfoil.data(:,1),xfoil.data(:,2), '--b', 'LineWidth', 2)
  legs_12 = [legs_12 'xfoil'];
  legend(legs_12, 'Interpreter', 'tex', 'FontSize', lfs)
end
 





