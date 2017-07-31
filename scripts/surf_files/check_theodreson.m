% Build a theoderson unsteady model

clear
clc
close all

% Load Xfoil data
re750 = importdata('test_ed36f128+14_re7.5e5.dat');


% Build parametric functions.
k0=0;
kmax=1000;
nk=kmax*100;
k=linspace(k0,kmax,nk);
k(1)=[];
C = besselh(1,2,k)./(besselh(1,2,k) + i*besselh(0,2,k));
F=real(C);
G=imag(C);

U0=1;
chord=1;
b=chord/2;
omega_k = k*U0/b;
a=-0.15;
ptch_amp=1.0*pi/180;
r2 = sqrt(1 + (a*k).^2);
added_m_eff_amp = -pi*ptch_amp*k.*r2;
added_m_eff_phi = asin(1./r2);

% This is the plot in Theoderson's paper.
figure(10)
plot(1./k,F); hold on; 
plot(1./k,-G, 'r'); 
grid on

figure(11)
plot(k,abs(added_m_eff_amp))
xlim([0 2])
ylabel('Added mass')

figure(12)
plot(k,added_m_eff_phi*180/pi)
xlim([0 2])
ylabel('Added mass phase shift')


% Building quasi-steady term
omegat = linspace(0,2*pi,100);
k=0.32;
U0=1.0;
chord=1.0;
b=chord/2;
omega=k*U0/b;
time=omegat/omega;
ptch_amp = 1.0*pi/180;               % radians
alpha0 = 3.4*pi/180;
a = -0.15;                    % pitch axis from mid chord

zi=sqrt(-1);
alpha = alpha0 + ptch_amp*exp(zi*omegat);
dalphadt = ptch_amp*zi*k*exp(zi*omegat);
lagterm = dalphadt*(1/2-a);
alpha_eff = alpha + lagterm;
Ck = besselh(1,2,k)./(besselh(1,2,k) + i*besselh(0,2,k));
alpha_eff_wake = real(alpha_eff*Ck);
r = sqrt(1 + 0.5-a);
ptch_amp_eff = ptch_amp*r*180/pi;

%% Added mass
r2 = sqrt(1 + (a*omega)^2);
m_eff = pi*ptch_amp*k*(cos(omegat) - a*k*sin(omegat));
m_eff_amp = -pi*ptch_amp*omega*r2;


% Total
cl_750k1 = interp1(re750.data(:,1),re750.data(:,2),real(alpha_eff_wake*180/pi));
cl_750k2r = interp1(re750.data(:,1),re750.data(:,2),real(alpha_eff*180/pi));
cl_750k2i = interp1(re750.data(:,1),re750.data(:,2),imag(alpha_eff*180/pi));

unsteady_cl_ideal = 2*pi*alpha_eff_wake + m_eff;
unsteady_cl_saab1 = cl_750k1 + m_eff;
unsteady_cl_saab2 = real(cl_750k2r) + m_eff;

alphareal = real(alpha);

figure(2)
plot(omegat,real(alpha_eff*180/pi)); hold on
plot(omegat,real(alpha_eff_wake*180/pi), 'r')
ylabel('\alpha_{eff}^{o}')
xlabel('\omega t')
grid on

figure(3)
plot(alphareal*180/pi,real(alpha_eff*180/pi)); hold on
plot(alphareal*180/pi,real(alpha_eff_wake*180/pi), 'r')
ylabel('\alpha_{eff}^{o}')
xlabel('\alpha^{o}')

figure(4)
plot(alphareal*180/pi,m_eff)
ylabel('M_{eff}')
xlabel('\alpha^{o}')

figure(5)
plot(alphareal*180/pi,unsteady_cl_saab1); hold on
plot(alphareal*180/pi,unsteady_cl_saab2, 'r');
ylabel('C_{l}')
xlabel('\alpha^{o}')


