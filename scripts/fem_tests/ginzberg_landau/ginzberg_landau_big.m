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
sigma = 0.05;
xofst = 20;
%u0 = normpdf(xgll-xofst,nu,sigma);
%u0 = 0.01*u0/max(max(u0));
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

ulag1 = 0*un;
ulag2 = 0*un;

mu_lag1 = 0*un;
mu_lag2 = 0*un;

mu_conv_lag1 = 0*un;
mu_conv_lag2 = 0*un;

mu_abs_lag1 = 0*un;
mu_abs_lag2 = 0*un;

cubic_lag1 = 0*un;
cubic_lag2 = 0*un;

quintic_lag1 = 0*un;
quintic_lag2 = 0*un;

ualag1 = 0*un;
ualag2 = 0*un;

conv_s = U0;
nek_linear_stiff = -conv_s*nek_conv - nek_lp;
%nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_linear_stiff;

% Parameters
deltat = 0.001;
istep = 0;
%nsteps = 500000;
time = 0;
iostep = 5000;
isave  = 500;
OMEGA=0.10;
Tosc=2*pi/OMEGA;
nsteps=ceil(4*Tosc/deltat);
verbose=1;
verbosestep=100;
ifpcg=0;

u_save = [];
t_save = []; 


%% SOLVE
h3=figure;
ax1=axes;
ax2=axes('XAxisLocation','Top');
axes(ax1)

nperiods=0;
for i = 1:nsteps
  t_it=tic;
  if (istep>0 && verbose && (mod(istep,verbosestep)==0))
    if (ifpcg)    
      disp(['Step, Time, Relative residual,iter,timing: ', num2str(istep), ', ' num2str(time), ', ' ...
        num2str(relres) ', ' num2str(pcgiter) ', ' num2str(time_it)]);
    else
      disp(['Step, Time, timing: ', num2str(istep), ', ' num2str(time), ', ' ...
        num2str(time_it)]);
    end
  end

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
  big_b=zeros(dof,1);
  gl_mass=zeros(dof,dof);
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

%   Constant Source term
    mu_term(:,:,els) = nek_intgd(:,:,els)*nek_mud(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
%   Time varying spatially constant source term
    mu_time = mut_0*sin(OMEGA*time);  
    mu_term(:,:,els) = mu_term(:,:,els) + nek_intgd(:,:,els)*mu_time*nek_intpm1d(:,:,els)*un(:,els);
    b31 = extk(2)*mu_term(:,els);
    b32 = extk(3)*mu_lag1(:,els);
    b33 = extk(4)*mu_lag2(:,els);
    b3 = deltat*(b31 + b32 + b33);

%   Time varying convectively unstable region
    mu_time = mut_conv*sin(OMEGA*time); 
    mu_conv_term(:,:,els) = mu_time*nek_intgd(:,:,els)*nek_mud_conv(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
    b41 = extk(2)*mu_conv_term(:,els);
    b42 = extk(3)*mu_conv_lag1(:,els);
    b43 = extk(4)*mu_conv_lag2(:,els);
    b4 = deltat*(b41 + b42 + b43);

%   Time varying absolute instability term
    mu_time = max([0 mut_abs*sin(OMEGA*time)]);
    mu_abs_term(:,:,els) = mu_time*nek_intgd(:,:,els)*nek_mutd(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
    b51 = extk(2)*mu_abs_term(:,els);
    b52 = extk(3)*mu_abs_lag1(:,els);
    b53 = extk(4)*mu_abs_lag2(:,els);
    b5 = deltat*(b51 + b52 + b53);

%   Cubic term
    intpd2_u = nek_intpm1d2(:,:,els)*un(:,els);
    cubic_term(:,:,els) = nek_intgd2(:,:,els)*(intpd2_u.^3); 
    b_cubic1 = extk(2)*cubic_term(:,els);
    b_cubic2 = extk(3)*cubic_lag1(:,els);
    b_cubic3 = extk(4)*cubic_lag2(:,els);
    b_cubic = deltat*(b_cubic1 + b_cubic2 + b_cubic3);

%   Quintic term
    quintic_intp = nek_intpm1d3(:,:,els)*un(:,els);
    quintic_term(:,:,els) = nek_intgd3(:,:,els)*(quintic_intp.^5); 
    b_quintic1 = extk(2)*quintic_term(:,els);
    b_quintic2 = extk(3)*quintic_lag1(:,els);
    b_quintic3 = extk(4)*quintic_lag2(:,els);
    b_quintic = deltat*(b_quintic1 + b_quintic2 + b_quintic3);

%   sum
    b(:,els) = -b1 + b2 + (b3 - b4 + b5) + b_cubic - b_quintic;

    gl_pos_j1 = (els-1)*N + 1;
    gl_pos_i1 = (els-1)*N + 1;

    gl_pos_j2 = (els)*N + 1;
    gl_pos_i2 = (els)*N + 1;
    
    gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = mass(:,:,els) + gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
    big_b(gl_pos_j1:gl_pos_j2) = b(:,els) + big_b(gl_pos_j1:gl_pos_j2);  

  end

%% Lag terms  
  ulag2 = ulag1;
  ulag1 = un;

  mu_lag2 = mu_lag1;
  mu_lag1 = mu_term;

  mu_conv_lag2 = mu_conv_lag1;
  mu_conv_lag1 = mu_conv_term;

  mu_abs_lag2 = mu_abs_lag1;
  mu_abs_lag1 = mu_abs_term;

  cubic_lag2 = cubic_lag1;
  cubic_lag1 = cubic_term; 

  quintic_lag2 = quintic_lag1;
  quintic_lag1 = quintic_term;
%%

  periodic = 0; 
  b = DSSUM(b,nels,periodic);

  mass(1,1) = 1;                % strong bc
  b(1,1)    = 2e-13 + 1e-13*rand(1);    % random noise O(1e-6)
%  b(1,1)    = 0;
%  b(1,1)    = 0.1*sin(time);
  big_b(1)  = b(1,1); 
%  mass(end,end) = 1;                % strong bc
%  b(end,end)    = 0;                % random noise O(1e-6) 

  gl_mass(1,:) = 0;
  gl_mass(1,1) = 1;
%  big_b(1) = 0;

  M=[];                         % preconditioner
  restrt=N*nels;
  max_it = 5;
  tol = 1e-9;

%  [soln err iter] = me_solve(mass,un,b,max_it,tol);
%  [soln err iter] = solve_strong(mass,un,b,max_it,tol,periodic);

  if ifpcg
    tol=1e-8;
    max_it=500;
    [v cflag relres pcgiter] = pcg(gl_mass,big_b,tol,max_it);
  else
    v = gl_mass\big_b;
  end  

  soln = gs_scatter(v,N,nels,periodic);

  un = soln;

  %% Just plot Points
  if mod(istep,iostep)==0    
    figure(h3)
    plot(ax1,xgll,un, '-k', 'MarkerSize', 16);
%    set(ax1,'Color', 'none');
%    xlim(ax1,[1 600])
    legend(ax1,['T/T_{osc}=' num2str(time/Tosc) '; istep=' num2str(istep) '; \mu_{t}=' num2str(mu_time)], 'Location', 'Best')
%    plot(ax2,mu_x,un, '-k', 'MarkerSize', 16);

    plot(ax2,xgll,mu_x+mut_0*sin(OMEGA*time)+max([0 sin(OMEGA*time)])*mu_xt, '--r', 'MarkerSize', 16);
    set(ax2,'XAxisLocation', 'Top')
    set(ax2,'YAxisLocation', 'Right')
%    set(ax2,'XDir', 'reverse')  
    set(ax2,'Color', 'none')
    ylabel(ax1,'A')
    ylabel(ax2,'\mu') 
    axes(ax2)
    pause(0.01)

  end

  if (mod(istep,isave)==0)
    u_save = [u_save un(:)];
    t_save = [t_save time];
  end       

  %dbstop in me_periodic at 335
  iperiod=floor(time/Tosc);    
  if (iperiod>nperiods)
    figure(10)
    surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
    view(2)
    pause(2)
    colorbar  
    nperiods=iperiod;  
  end      
         
  time_it = toc(t_it);
end

figure(10)
surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
view(2)

save('homogeneous_abs1.mat')

