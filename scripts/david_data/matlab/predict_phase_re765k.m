% Trying to predict phase plots from previous data

clear
clc
close all

run get_predictions

onlymean = 0;

static_model = load('14_static_models_765k.mat');
steady_shift = 0.;
uoo=24;

%% ideal numbers
mean_aoa=3.2;
dalpha=1.0;
k=0.4;           % min value is 0.02
%%

%% closest to data
%mean_aoa=3.2;
%dalpha=0.94;
%k=0.32;           % min value is 0.02

%% propose
mean_aoa=3.4;
dalpha=1.0;
k=0.4;           % min value is 0.02


U0=1.;
c=1.0;
semichord=c/2;
omega=k*U0/semichord;

omegat = linspace(0,4*pi,1000);
time=omegat/omega;
intg_const = interp1(movk,movintg,k,'linear');
intg_const = intg_const*dalpha*pi/180;
intg_min = intg_const*0.85;
intg_max = intg_const*1.15;

norm_intg_const = interp1(movk,movintgnorm,k,'linear');
intg_const = norm_intg_const*dalpha*pi/180*30/uoo;
intg_min = intg_const*0.9;
intg_max = intg_const*1.1;

philag = interp1(movk,movphi,k,'linear');
philag_min=0.8*philag;
philag_max=1.2*philag;

gammalag = interp1(movk,movgamma,k,'linear');
gammalag_min=0.7*gammalag;
gammalag_max=1.3*gammalag;

[intg_const philag*180/pi gammalag*180/pi]

aoa = mean_aoa + dalpha*sin(omegat);
pitch = intg_const*cos(omegat);
pitch_min = intg_min*cos(omegat);
pitch_max = intg_max*cos(omegat);

alphalagg =  mean_aoa + dalpha*sin(omegat + philag);
alphalagg_min =  mean_aoa + dalpha*sin(omegat + philag_min);
alphalagg_max =  mean_aoa + dalpha*sin(omegat + philag_max);

cz_lag = interp1(static_model.alpha,static_model.cz,alphalagg,'linear');
cz_lagmin = interp1(static_model.alpha,static_model.cz,alphalagg_min,'linear');
cz_lagmax = interp1(static_model.alpha,static_model.cz,alphalagg_max,'linear');

delta_static_cz = max(cz_lag)-min(cz_lag);
disp(['Max change in static Cz= ' num2str(delta_static_cz)])

cz = pitch + cz_lag;
cz_1 = pitch_min+cz_lagmin;
cz_2 = pitch_max+cz_lagmax;
cz_3 = pitch_min+cz_lagmax;
cz_4 = pitch_max+cz_lagmin;

col1 = lines(5);
if (onlymean)
  legs_12 = {'Cz' 'Static Cz'};
else
  legs_12 = {'Cz' 'Cz1' 'Cz2' 'Cz3' 'Cz4' 'Static Cz'};
end  

figure(20)
plot(aoa,cz, 'LineWidth', 4, 'Color', 'k'); hold on
if (~onlymean)
  plot(aoa,cz_1, '-.', 'LineWidth', 2, 'Color', col1(2,:));
  plot(aoa,cz_2, '-.', 'LineWidth', 2, 'Color', col1(3,:));
  plot(aoa,cz_3, '-.', 'LineWidth', 2, 'Color', col1(4,:));
  plot(aoa,cz_4, '-.', 'LineWidth', 2, 'Color', col1(5,:));
end 
plot(static_model.alpha,static_model.cz, '-b', 'LineWidth', 2)
legend(legs_12, 'Interpreter', 'tex', 'FontSize', lfs)
xmin = mean_aoa-dalpha*1.4;
xmax = mean_aoa+dalpha*1.4;
xlim([xmin xmax])
title('Phase Potrait')
ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
grid on
filename=['model_cz_alpha.eps'];
SaveFig(gcf,filename, destn, ifcols)

if (onlymean)
  legs_13 = {'Cz'};
else  
  legs_13 = {'Cz' 'Cz1' 'Cz2' 'Cz3' 'Cz4'};
end

figure(21)
plot(time,cz, 'LineWidth', 3, 'Color', 'k'); hold on
if (~onlymean)
  plot(time,cz_1, '-.', 'LineWidth', 2, 'Color', col1(2,:));
  plot(time,cz_2, '-.', 'LineWidth', 2, 'Color', col1(3,:));
  plot(time,cz_3, '-.', 'LineWidth', 2, 'Color', col1(4,:));
  plot(time,cz_4, '-.', 'LineWidth', 2, 'Color', col1(5,:));
end
legend(legs_13, 'Interpreter', 'tex', 'FontSize', lfs)
title('Cz-time')
ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('t', 'Interpreter', 'tex', 'FontSize', fs)
filename=['model_cz_time.eps'];
SaveFig(gcf,filename, destn, ifcols)


%% Xfoil
ifxfoil = 0;
if ifxfoil
  xfoil = importdata('polar_re765k_ed36f128+14.dat');
  figure(20)
  plot(xfoil.data(:,1),xfoil.data(:,2), '--b', 'LineWidth', 2)
  legs_12 = [legs_12 'xfoil'];
  legend(legs_12, 'Interpreter', 'tex', 'FontSize', lfs)
end
 





