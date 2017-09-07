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
r2 = sqrt(1 + (a*omega_k).^2);
added_m_eff_amp = -pi*ptch_amp*omega_k.*r2;


% This is the plot in Theoderson's paper.
figure(10)
plot(1./k,F); hold on; 
plot(1./k,-G, 'r'); 
grid on

figure(11)
plot(k,abs(added_m_eff_amp))
xlim([0 2])
ylabel('Added mass')




% Building quasi-steady term
omegat = linspace(0,2*pi,100);
k=0.4;
U0=1.0;
chord=1.0;
b=chord/2;
omega=k*U0/b;
time=omegat/omega;
ptch_amp = 1.0*pi/180;               % radians
alpha0 = 3.4*pi/180;
a = -0.15;                    % pitch axis from mid chord

alpha = alpha0 + ptch_amp*sin(omegat);
dalphadt = ptch_amp*omega*cos(omegat);
lagterm = dalphadt*(1/2-a);
alpha_eff = alpha + lagterm;
Ck = besselh(1,2,k)./(besselh(1,2,k) + i*besselh(0,2,k));
alpha_eff_wake = real(alpha_eff*Ck);
r = sqrt(1 + 0.5-a);
ptch_amp_eff = ptch_amp*r*180/pi;

%% Added mass
r2 = sqrt(1 + (a*omega)^2);
m_eff = pi*ptch_amp*omega*(cos(omegat) - a*omega*sin(omegat));
m_eff_amp = -pi*ptch_amp*omega*r2;


% Total
cl_750k = interp1(re750.data(:,1),re750.data(:,2),alpha_eff_wake*180/pi);
unsteady_cl_ideal = 2*pi*alpha_eff_wake + m_eff;
unsteady_cl_saab = cl_750k + m_eff;


figure(2)
plot(omegat,alpha_eff*180/pi)
ylabel('\alpha_{eff}^{o}')
xlabel('\omega t')
grid on

figure(3)
plot(alpha*180/pi,alpha_eff*180/pi)
ylabel('\alpha_{eff}^{o}')
xlabel('\alpha^{o}')

figure(4)
plot(alpha*180/pi,m_eff)
ylabel('M_{eff}')
xlabel('\alpha^{o}')

figure(5)
plot(alpha*180/pi,unsteady_cl_saab)
ylabel('C_{l}')
xlabel('\alpha^{o}')


