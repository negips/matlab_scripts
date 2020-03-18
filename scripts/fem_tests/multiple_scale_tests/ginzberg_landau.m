% Small test simulation.

clear
clc
close all

% Initialize 1D spectral element
spec_element_init

% Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];

who

x=xall;
x2 = linspace(d_start,d_end,interp_pts);

% use a normalized gaussian
nu = 0;
sigma = 0.2;
xofst = 50;
%u0 = normpdf(xgll-xofst,nu,sigma);
%u0 = 0.1*u0/max(max(u0));
%u0_int = normpdf(x2,mu,sigma);
%u0_int = u0_int/max(u0_int);

u0 = 0.000*sin(0.01*pi*xgll);
%u0_int = 0.000*sin(0.01*pi*x2);
%u0(50) = 1.;

un = u0;

%% testing
%h_isoln = figure;
%plot(x2,u0_int,'o');
%title('initial solution');

%u0_deriv = pi*0.1*cos(pi*x2);

deltat = 0.001;
istep = 0;
nsteps = 10000;

ulag1 = 0*un;
ulag2 = 0*un;

mu_lag1 = 0*un;
mu_lag2 = 0*un;

nl_lag1 = 0*un;
nl_lag2 = 0*un;

ualag1 = 0*un;
ualag2 = 0*un;

conv_s = 1/sqrt(3);
nek_linear_stiff = -conv_s*nek_conv - nek_lp;
%nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_linear_stiff;

time = 0;

%% SOLVE
h3=figure;

for i = 1:nsteps

  istep = istep+1;
  time = time+deltat; 

  % set time integration operators
  if istep == 1
       bdfk = bdf1;
       extk = ex0;
  elseif istep == 2
       bdfk = bdf2;
       extk = ex1;
  else
       bdfk = bdf3;
       extk = ex2;
  end

  b = 0*un;
  soln = 0*b;
  relres=zeros(nels);

  for els=1:nels
    mass(:,:,els) = nek_mass(:,:,els)*bdfk(1);

    b11 = bdfk(2)*un(:,els); 
    b12 = bdfk(3)*ulag1(:,els); 
    b13 = bdfk(4)*ulag2(:,els);
    b1 = nek_mass(:,:,els)*(b11+b12+b13);

%   Conv + dissipation
    stiff(:,:,els) = nek_linear_stiff(:,:,els);
    b21 = extk(2)*un(:,els);
    b22 = extk(3)*ulag1(:,els);
    b23 = extk(4)*ulag2(:,els);
    b2 = deltat*stiff(:,:,els)*(b21 + b22 + b23);

%   Source term
    mu_term(:,:,els) = nek_intgd(:,:,els)*nek_mud(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
    b31 = extk(2)*mu_term(:,els);
    b32 = extk(3)*mu_lag1(:,els);
    b33 = extk(4)*mu_lag2(:,els);
    b3 = deltat*(b31 + b32 + b33);

%   Cubic term
    intpd2_u = nek_intpm1d2(:,:,els)*un(:,els);
    nl_term(:,:,els) = nek_intgd2(:,:,els)*(intpd2_u.^3); 
    b41 = extk(2)*nl_term(:,els);
    b42 = extk(3)*nl_lag1(:,els);
    b43 = extk(4)*nl_lag2(:,els);
    b4 = deltat*(b41 + b42 + b43);

%   sum
    b(:,els) = -b1 + b2 + b3 - b4;

  end
  mu_lag2 = mu_lag1;
  mu_lag1 = mu_term;

  nl_lag2 = nl_lag1;
  nl_lag1 = nl_term; 

  periodic = 0; 
  b = DSSUM(b,nels,periodic);

  mass(1,1) = 1;                % strong bc
%  b(1,1)    = -1e-3 + 2e-3*rand(1);    % random noise O(1e-6)
  b(1,1)    = 0;
  b(1,1)    = 0.001*sin(2*pi*time); 
%  mass(end,end) = 1;                % strong bc
%  b(end,end)    = 0;                % random noise O(1e-6) 

  M=[];                         % preconditioner
  restrt=N*nels;
  max_it = 5;
  tol = 1e-9;

%  [soln err iter] = me_solve(mass,un,b,max_it,tol);
  [soln err iter] = solve_strong(mass,un,b,max_it,tol,periodic);

  ulag2 = ulag1;
  ulag1 = un;
  un = soln;

  %% Just plot Points
  figure(h3)
  plot(xgll,un, '-k', 'MarkerSize', 16);
  % hold on

  %dbstop in me_periodic at 335

end



