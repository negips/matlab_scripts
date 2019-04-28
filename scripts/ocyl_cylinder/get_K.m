%  Get Stiffness values for different reduced frequencies

m = 10;     % mass ratio
D = 1.0;    % Diameter
r = D/2;    % radius
inertia = m*pi*r^2;

omega_vortex_shedding = 0.7665486;

omega_ratio = [0.1 0.5 1 5 10];

w_struc = omega_ratio*omega_vortex_shedding;

stiffness = inertia*(w_struc.^2);
stiffness'

[omega_ratio' stiffness']





