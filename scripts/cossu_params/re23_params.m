% Evaluate parameters for the simulation

clear
clc

% Free-stream velocity
U0 = 1.0;
% Cylinder Diameter
D=1.0;

Re1=23.512;

% Cylinder area
A0 = 1/4*pi*D^2;
%A0 = 1;

%least_damped_omega = 0.0977219;
cossu_omegan = 1.326038;

ratio = 0.6;
%omegan = least_damped_omega*ratio;
omegan = cossu_omegan*ratio;

rho_star = 7.0;

% Density ratio (solid to fluid)
m_star = rho_star*A0;
K_star = m_star*omegan^2;
critical_damping = 2*m_star*omegan;
damping_ratio = 0.00754126;

D_star = damping_ratio*critical_damping;

% Spring-mass frequency
fn = omegan/(2*pi);

% Reduced Velocity
U_star = (U0/D)*(1/fn);
omegad = omegan*sqrt(1-damping_ratio^2);
lambda_s = -omegan*damping_ratio + 1i*omegad;

disp(' ')

disp(['Ratio:                             ', num2str(ratio,4)])
disp(['Density Ratio:                     ', num2str(rho_star,7)])
disp(['Mass:                              ', num2str(m_star,7)])
disp(['Undamped Angular Frequency:        ', num2str(omegan,7)])
disp(['Structural Angular Frequency:      ', num2str(omegad,7)])
disp(['Reduced Frequency:                 ', num2str(U_star,7)])
disp(['Normalized stiffness:              ', num2str(K_star,7)])
disp(['Normalized Damping:                ', num2str(D_star,7)])
disp(['Critical Damping coefficient:      ', num2str(critical_damping,7)])
disp(['Damping Ratio:                     ', num2str(damping_ratio,7)])
disp(['Structural Eigenvalue:             ', num2str(real(lambda_s),7) ' + i',num2str(imag(lambda_s),7) ])






