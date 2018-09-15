% Random cl response

clear
clc
close all

phi = linspace(0,6*pi,1000);
dalpha = 1.0*pi/180;
phi_lag = 40*pi/180;
phi_gain = 120*pi/180;
amp = 1.0;
nl_amp = 0.5;

a0 = 10*pi/180;

a  = a0+dalpha*sin(phi);
a1 = a0+dalpha*sin(phi+phi_lag);
a2 = a0+dalpha*sin(phi+phi_gain);
nl = dalpha*2*sin(2*(phi+phi_lag));

cl_nl = pi*nl;
cl = 2*pi*(a1) - nl_amp*(nl) + amp*(a2);

plot(a,cl)
