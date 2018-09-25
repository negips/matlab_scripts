% just testing things right now

clear
clc
close all

alpha_deg = linspace(-18,18,200);
alpha = alpha_deg*pi/180;

rt1 = 0.;         % 0 degrees
rt2 = 17*pi/180;  % 2nd root in degrees
rt3 = -rt2;
grad_0 = -2*pi/(rt1*rt2 + rt1*rt3 + rt2*rt3);         % gradient at alpha=0;


cl = grad_0*(rt1 - alpha).*(rt2 - alpha).*(rt3 - alpha);

plot(alpha_deg,cl)
grid on
