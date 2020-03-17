% Evaluate parameters for the simulation

clear
clc

% Free-stream velocity
U0 = 1.0;

Re=60.0;

rho_star = 10.0;


% Ellipse axes (radii)
a=0.25;           % x-radius
b=0.5;            % y-radius

D = 2*b;

% Ellipse area
A0 = pi*(a*b);
M0 = rho_star*A0;
I0 = M0/4*(a^2 + b^2);

Ustar  = 7.0;
Fn     = 1/Ustar;
fn     = Fn*U0/D;
omegan = 2*pi*fn;

% Density ratio (solid to fluid)
m_star = I0;
K_star = m_star*omegan^2;
critical_damping = 2*m_star*omegan;
damping_ratio = 0.0;

D_star = damping_ratio*critical_damping;

% Spring-mass frequency
fn = omegan/(2*pi);

% Reduced Velocity
omegad = omegan*sqrt(1-damping_ratio^2);
lambda_s = -omegan*damping_ratio + 1i*omegad;

disp(' ')

disp(['Density Ratio:                     ', num2str(rho_star,7)])
disp(['Inertia:                           ', num2str(m_star,7)])
disp(['Undamped Angular Frequency:        ', num2str(omegan,7)])
disp(['Structural Angular Frequency:      ', num2str(omegad,7)])
disp(['Reduced Frequency:                 ', num2str(Ustar,7)])
disp(['Normalized stiffness:              ', num2str(K_star,7)])
disp(['Normalized Damping:                ', num2str(D_star,7)])
disp(['Critical Damping coefficient:      ', num2str(critical_damping,7)])
disp(['Damping Ratio:                     ', num2str(damping_ratio,7)])
disp(['Structural Eigenvalue:             ', num2str(real(lambda_s),7) ' + i',num2str(imag(lambda_s),7) ])






