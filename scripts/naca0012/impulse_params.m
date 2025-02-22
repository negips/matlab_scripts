% Parameters for the Poirel (2008) Self-sustained aeroelastic oscillations of a NACA0012 
%     airfoil at low-to-moderate Reynolds numbers

clear
clc

x_theta = 0.15;               % Distance between Elastic axis and centre of mass
a_h     = -0.63;              % Non-dimensional location of elastic axis from the mid-chord (Non-dimensionalized by half-chord)
I_s     = 0.00135;            % Mass moment of Inertial of rotating parts (Kg-m^2)
m_theta = 0.771;              % Mass of parts in rotation (Kg)
K_s     = 0.30;               % Structural Stiffness (N-m/rad)
D_s     = 0.002;              % Approx Struc. Dissipation (N-m-s)
c       = 0.156;              % Chord (m)
S       = 0.61;               % Span  (m)


U_inf   = 7.5;                % Free-stream velocity (m/s)
Re      = 7.7E+04;            % Reynolds Number
nu      = U_inf*c/Re;         % Kinematic viscosity (m^2/s)


air = importdata('air_prop.dat');

rho = interp1(air.data(:,8),air.data(:,9),nu);        % Density (kg/m^3)
mu  = interp1(air.data(:,8),air.data(:,5),nu);        % Dynamic viscosity (kg/(m-s))
T   = interp1(air.data(:,8),air.data(:,1),nu);        % Temperature (K)

nu = mu/rho;
time_scale = 1.0;
U_scale = c/time_scale;

Re = U_scale*c/nu;

axis_x0 = c/2 + a_h*c/2;
X0_star = axis_x0/c;

% Non-dimensionalization (Rotational)

I_star = I_s/(rho*c^5);

D_star = D_s/(rho*U_scale*c^4);

K_star = K_s/(rho*U_scale^2*c^3);

disp(['Reynolds No (Re) :          ', num2str(Re,7)])
disp(['Temperature (T)  :          ', num2str(T,7)])
disp(['Density     (rho):          ', num2str(rho,7)])
disp(['Inertia     (I*) :          ', num2str(I_star,7)])
disp(['Struc Dissip(D*) :          ', num2str(-D_star,7)])
disp(['Stiffness   (K*) :          ', num2str(K_star,7)])
disp(['Axis       (X0*) :          ', num2str(X0_star,7)])














