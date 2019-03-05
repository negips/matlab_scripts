% Small test simulation.

clear
clc
close all

% Initialize 1D spectral element
conv_diff_init

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

u0 = 0.1*sin(2*pi*xgll);
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

mu_t_lag1 = 0*un;
mu_t_lag2 = 0*un;

nl_lag1 = 0*un;
nl_lag2 = 0*un;

quintic_lag1 = 0*un;
quintic_lag2 = 0*un;

ualag1 = 0*un;
ualag2 = 0*un;

conv_s = 1;
Re=100;
nek_linear_stiff = -conv_s*nek_conv - 0*nek_lp;
nek_linear_stiff_full  = -conv_s*nek_convd;
%nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_linear_stiff;

% Parameters
deltat = 0.01;
istep = 0;
%nsteps = 500000;
time = 0;
iostep = 100;
isave  = 200;
OMEGA=1.0;
Tosc=2*pi/OMEGA;
amp=mut;
nsteps=ceil(500*Tosc/deltat);

u_save = [];
t_save = []; 

%% SOLVE
h3=figure;
ax1=axes;
%ax2=axes('XAxisLocation','Top');
%axes(ax1)


%% Statistics
stat_alpha=0;
stat_beta=0;
stat_atime=0;
umean = 0*un;
u2mean = 0*un;

nperiods=0;
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
  big_b=zeros(dof,1);
  gl_mass=zeros(dof,dof);
  relres=zeros(nels);

  for els=1:nels
    mass(:,:,els) = nek_mass(:,:,els)*bdfk(1);
%    mass(:,:,els) = nek_intgd(:,:,els)*nek_intpm1d(:,:,els)*bdfk(1);  

    b11 = bdfk(2)*un(:,els); 
    b12 = bdfk(3)*ulag1(:,els); 
    b13 = bdfk(4)*ulag2(:,els);
    b1 = nek_mass(:,:,els)*(b11+b12+b13);
%    b1  = nek_intgd(:,:,els)*nek_intpm1d(:,:,els)*(b11+b12+b13); 

%   Conv + dissipation
    stiff(:,:,els) = nek_linear_stiff(:,:,els);
%    stiff(:,:,els) = nek_linear_stiff_full(:,:,els);  
    b21 = extk(2)*un(:,els);
    b22 = extk(3)*ulag1(:,els);
    b23 = extk(4)*ulag2(:,els);
    b2 = deltat*stiff(:,:,els)*(b21 + b22 + b23);

%   Constant Source term
%    mu_term(:,:,els) = nek_intgd(:,:,els)*nek_mud(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
%    b31 = extk(2)*mu_term(:,els);
%    b32 = extk(3)*mu_lag1(:,els);
%    b33 = extk(4)*mu_lag2(:,els);
%    b3 = deltat*(b31 + b32 + b33);

%   Cubic term
%    intpd2_u = nek_intpm1d2(:,:,els)*un(:,els);
%    nl_term(:,:,els) = nek_intgd2(:,:,els)*(intpd2_u.^3); 
%    b41 = extk(2)*nl_term(:,els);
%    b42 = extk(3)*nl_lag1(:,els);
%    b43 = extk(4)*nl_lag2(:,els);
%    b4 = deltat*(b41 + b42 + b43);

%   time modulated Source term
%    mu_time = amp*sin(OMEGA*time);
%    mu_t_term(:,:,els) = mu_time*nek_intgd(:,:,els)*nek_mutd(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
%    b51 = extk(2)*mu_t_term(:,els);
%    b52 = extk(3)*mu_t_lag1(:,els);
%    b53 = extk(4)*mu_t_lag2(:,els);
%    b5 = deltat*(b51 + b52 + b53);

%   Quintic term
%    quintic_intp = nek_intpm1d3(:,:,els)*un(:,els);
%    quintic_term(:,:,els) = nek_intgd3(:,:,els)*(quintic_intp.^5); 
%    b61 = extk(2)*quintic_term(:,els);
%    b62 = extk(3)*quintic_lag1(:,els);
%    b63 = extk(4)*quintic_lag2(:,els);
%    b6 = deltat*(b61 + b62 + b63);

%   sum
%    b(:,els) = -b1 + b2 + b3 + b4 + b5 - b6;
    b(:,els) = -b1 + b2;  

    gl_pos_j1 = (els-1)*N + 1;
    gl_pos_i1 = (els-1)*N + 1;

    gl_pos_j2 = (els)*N + 1;
    gl_pos_i2 = (els)*N + 1;
    
    gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = mass(:,:,els) + gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
    big_b(gl_pos_j1:gl_pos_j2) = b(:,els) + big_b(gl_pos_j1:gl_pos_j2);  

  end

  ulag2 = ulag1;
  ulag1 = un;

%  mu_lag2 = mu_lag1;
%  mu_lag1 = mu_term;
%
%  mu_t_lag2 = mu_t_lag1;
%  mu_t_lag1 = mu_t_term;
%
%  nl_lag2 = nl_lag1;
%  nl_lag1 = nl_term; 
%
%  quintic_lag2 = quintic_lag1;
%  quintic_lag1 = quintic_term; 

  periodic = 1; 
  b = DSSUM(b,nels,periodic);

  if ~periodic
    mass(1,1) = 1;                % strong bc
%    b(1,1)    = 2e-4 + 1e-4*rand(1);    % random noise O(1e-6)
%    b(1,1)    = 0;
    b(1,1)    = 0.1*sin(time);
    big_b(1)  = b(1,1); 
%    mass(end,end) = 1;                % strong bc
%    b(end,end)    = 0;                % random noise O(1e-6) 

    gl_mass(1,:) = 0;
    gl_mass(1,1) = 1;
%    big_b(1) = 0;
  else
    big_b(1) = big_b(1)+big_b(end);
    big_b(end)= [];  
    gl_mass_temp = gl_mass;  
    gl_mass(1,:)=gl_mass(1,:)+gl_mass_temp(end,:);
    gl_mass(:,1)=gl_mass(:,1)+gl_mass(:,end);
    gl_mass(end,:)=[];
    gl_mass(:,end)=[];    
  end              

  M=[];                         % preconditioner
  restrt=N*nels;
  max_it = 5;
  tol = 1e-9;

%  [soln err iter] = me_solve(mass,un,b,max_it,tol);
%  [soln err iter] = solve_strong(mass,un,b,max_it,tol,periodic);

  v = gl_mass\big_b;
  v(end+1)=v(1);    

  soln = gs_scatter(v,N,nels,periodic);

  un = soln;

  %% Just plot Points
  if mod(istep,iostep)==0    
    figure(h3)
    plot(ax1,xgll,un, '-k', 'MarkerSize', 16);
%    xlim(ax1,[1 500])
    legend(['T=' num2str(time) '; istep=' num2str(istep)], 'Location', 'Best')
    ylim([-0.1 0.1])  
%    plot(ax2,mu_x,un, '-k', 'MarkerSize', 16);
%    set(ax2,'XAxisLocation', 'Top')
%    xlim(ax2,[min(mu_x(:)) max(mu_x(:))]+mu_time)  
%    xticks2 = get(ax2,'XTick');
%    set(ax2,'XDir', 'reverse')  
%    set(ax2,'XTick',xticks2(end:-1:1))
%    set(ax2,'Color', 'none')
%    set(ax2,'YTick',[])
    ylabel(ax1,'A')    
    pause(0.001)

  end

  if (mod(istep,isave)==0)
    u_save = [u_save un(:)];
    t_save = [t_save time];
  end       

%  iperiod=floor(time/Tosc);    
%  if (iperiod>nperiods)
%    figure(10)
%    surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
%    view(2)
%    pause(2)
%    colorbar  
%    nperiods=iperiod;  
%  end      

  stat_atime=stat_atime+deltat;
  stat_beta=deltat/stat_atime;
  stat_alpha = 1 - stat_beta;

  umean = stat_alpha*umean + stat_beta*un;
  u2mean = stat_alpha*u2mean + stat_beta*(un.^2);

end

%figure(10)
%surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
%view(2)

figure(2)
urms=u2mean-umean.^2;
plot(xgll,urms)

save('conv_diag2.mat')


