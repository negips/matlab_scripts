% Calculate simulation parameters for the K. Menon, R. Mittal. "Computational Modeling and Analysis of Aeroelastic Wing Flutter", 2018 Fluid Dynamics Conference, AIAA Aviation forum (AIAA 2018-3080)

clear
clc
close all

Uinf = 1.0;
c = 1.0;
rho = 1.0;

U_star = 7.5;


I_n = 2.07;

I = I_n*rho*c^4/2;
fn = Uinf/c/U_star;
k_theta = (2*pi*fn)^2*I;
k_theta_n = 2*k_theta/(rho*Uinf^2*c^2);
b_theta_n = 0.3*sqrt(I_n*k_theta_n);

disp(['Normalized Inertia:     ',num2str(I_n,5)])
disp(['Normalized Stiffness:   ',num2str(k_theta_n,5)])
disp(['Normalized Damping:     ',num2str(b_theta_n,5)])

