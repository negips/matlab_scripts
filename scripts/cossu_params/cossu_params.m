% Evaluate parameters for the simulation

clear
clc

% Free-stream velocity
U0 = 1.0;
% Cylinder Diameter
D=1.0;

Re1=23.512;
Re2=47.024;

Re=Re1;
caseno=1;

tol=1.0e-7;
if abs(Re-Re1)<tol
  switch (caseno)
    case 1
      Lambda_s = -5.0e-3 + 1i*0.663;
      n = 1/70;
    case 2
      Lambda_s = -5.0e-3 + 1i*0.663;
      n = 1/7;
  end

else abs(Re-Re2)<tol

  switch (caseno)
    case 1
      Lambda_s = -5.0e-3 + 1i*0.663;
      n=1./7000;
    case 2
      Lambda_s = -5.0e-3 + 1i*0.663;
      n=1./700;
  end
end

% Cylinder area
A0 = 1/4*pi*D^2;
%A0 = 1;

% Density ratio (solid to fluid)
m_star = 1/n*A0;
Gamma = -2*real(Lambda_s);
damping = m_star*Gamma;

omega_damped = imag(Lambda_s);
omega_damped_squared = omega_damped^2;

omega_undamped_squared = omega_damped_squared + (Gamma/2)^2;
omega_undamped = sqrt(omega_undamped_squared);

K_star = m_star*omega_undamped_squared;           % Stiffness

critical_damping = 2*sqrt(K_star*m_star);         % 2*sqrt(k*m) Critical Damping coefficient

%Kai_star = Gamma/2/omega_undamped;              % Damping
Kai_star = Gamma*m_star;              % Damping

% Spring-mass frequency
fn = omega_undamped/(2*pi);

% Reduced Velocity
U_star = (U0/D)*(1/fn);

disp('Cossu & Morino (2000) On the instability of a spring-mounted circular cylinder in a viscous flow at low Reynolds Numbers. J. Fluids and Struc.')
disp('  ')

disp(['Density Ratio:                     ', num2str(1/n,7)])
disp(['Mass:                              ', num2str(m_star,7)])
disp(['Undamped Angular Frequency:        ', num2str(omega_undamped,7)])
disp(['Structural Angular Frequency:      ', num2str(omega_damped,7)])
disp(['Reduced Frequency:                 ', num2str(U_star,7)])
disp(['Normalized stiffness:              ', num2str(K_star,7)])
disp(['Normalized Damping:                ', num2str(Kai_star,7)])
disp(['Critical Damping coefficient:      ', num2str(critical_damping,7)])

lambda_s = -1/2*(Kai_star/m_star) + sqrt(-1)*sqrt(K_star/m_star - 1/4*(Kai_star/m_star)^2 )






