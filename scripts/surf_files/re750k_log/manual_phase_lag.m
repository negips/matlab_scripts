% Check Phase lag manually

clear
clc
close all

xfoil = importdata('test_ed36f128+14_re7.5e5.dat');
xfoilalpha = xfoil.data(:,1);
xfoilcl = xfoil.data(:,2);
figure(2)
plot(xfoilalpha,xfoilcl, '--k')
hold on


phi=-100*pi/180;
intg_const=0.15;
phase_lag = -pi/2;
dalpha = 1.0;
alpha0 = 3.4;
theta = -50*pi/180; 

omegat = linspace(0,2*pi,100);    

alpha = alpha0 + dalpha*sin(omegat + phase_lag);

pitch = dalpha*sin(omegat + phi + phase_lag);
alpha_lagg = alpha0 + pitch;

added_mass = intg_const*cos(omegat + theta + phase_lag);

cl_lagg = interp1(xfoilalpha,xfoilcl,alpha_lagg,'linear');

cl_pred = 0 + added_mass + cl_lagg;

figure(2)
plot(alpha,cl_pred)    

