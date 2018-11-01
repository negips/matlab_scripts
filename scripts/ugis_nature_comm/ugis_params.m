% Evaluate parameters for the simulation

clear
clc

% Free-stream velocity
U0 = 1.0;

% Cylinder Diameter
D=1.0;

% Reduced velocity
U_star = inf;

% Dimensionless Damping ratio
Xi_star = 0.000;

% Density ratio
m_star = 1.01;

% Splitter plate dimensions
a = D;            % Length
b = 0.02*D;       % Thickness

% Area moment of inertia (Splitter plate) about splitter plate centre of mass
A_sp = a*b;                   % Area of splitter plate
I_sp = 1/12*A_sp*(a^2 + b^2);

% Area moment about centre of cylinder
dr = 1.0;         % distance of centre of mass of plate from (0,0)
AI_sp_00 = I_sp + A_sp*dr^2;

% Area moment of inertia of cylinder about (0,0)
A_cyl = pi*(D/2)^2;
AI_cyl_00 = pi/2*(D/2)^4;

% Total Area moment of inertia about (0,0)
AI_00 = AI_cyl_00 + AI_sp_00;

% Mass moment of inertia about (0,0)
MI_00 = m_star*AI_00;

% Simulation frequency scale
Ff = U0/D;

% Reduced frequency
Fn = 1/U_star;          % dimensionless

% Spring-mass frequency
fn = Fn*(U0/D);         % dimensional

% Non-dimensional stiffness
K_star = MI_00*(2*pi*Fn)^2;

% Dimensional Damping ratio
Kai_star = Xi_star*2*sqrt(K_star*MI_00);

disp('Lacis et al. (2014) Passive appendeges generate drift through symmetry breaking. Nature Communications')
disp('  ')

disp(['Density Ration:              ', num2str(m_star,7)])
disp(['Reduced Frequency:           ', num2str(U_star,7)])
disp(['Angular Frequency:           ', num2str(sqrt(K_star/MI_00),7)])
disp(['Mass moment of inertia:      ', num2str(MI_00,7)])
disp(['Normalized stiffness:        ', num2str(K_star,7)])
disp(['Damping ratio:               ', num2str(Xi_star,7)])
disp(['Normalized Damping ratio:    ', num2str(Kai_star,7)])








