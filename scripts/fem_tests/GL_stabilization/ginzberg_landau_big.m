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

u0 = rand([length(xall),1]);
% use a normalized gaussian
nu = 0;
sigma = 0.2;
xofst = 5;
u0 = 2*normpdf(xgll-xofst,nu,sigma);
u0 = 1e-6*u0/max(max(u0));

un = u0;

%% testing
%h_isoln = figure;
%plot(x2,u0_int,'o');
%title('initial solution');

%u0_deriv = pi*0.1*cos(pi*x2);

ulag1 = 0*un;
ulag2 = 0*un;

source_lag1 = 0*un;
source_lag2 = 0*un;

cubic_lag1 = 0*un;
cubic_lag2 = 0*un;

quintic_lag1 = 0*un;
quintic_lag2 = 0*un;

forcing_lag1 = 0*un;
forcing_lag2 = 0*un;

conv_s = U0;
nek_linear_stiff = -conv_s*nek_conv - nek_lp;
%nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_linear_stiff;

% Parameters
deltat = 0.002;
istep = 0;
%nsteps = 500000;
time = 0;
iostep = 10;
isave  = 100;
OMEGA=0.01;
Tosc=2*pi/OMEGA;
nsteps=ceil(1.5*Tosc/deltat);
surfupd = 1.00;               % make surface evry surfupd of a period
verbose=1;
verbosestep=500;
ifpcg=0;
ifdealias=0;
ifrenorm=0;

u_save = [];
u_renorm = [];
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
    if ifdealias  
      mu_term(:,els) = nek_intgd(:,:,els)*nek_mud(:,:,els)*nek_intpm1d(:,:,els)*un(:,els);
    else  
      mu_term(:,els) = nek_mass(:,:,els)*nek_mu(:,:,els)*un(:,els);
    end  

%   Convectively unstable region (constant in time)
    mu_c = mu1;
    if ifdealias 
      mu_conv_term(:,els) = mu_c*nek_intgd(:,:,els)*nek_mud_conv(:,:,els)*nek_intpm1d(:,:,els)*un(:,els); 
    else  
      mu_conv_term(:,els) = mu_c*nek_mass(:,:,els)*nek_mu_conv(:,:,els)*un(:,els); 
    end

    source_term(:,els) = mu_term(:,els); 
    b_source = deltat*(extk(2)*source_term(:,els) + extk(3)*source_lag1(:,els) + extk(4)*source_lag2(:,els));


%   Cubic term
    if ifdealias  
      intpd2_u = nek_intpm1d2(:,:,els)*un(:,els);
      cubic_term(:,els) = nek_intgd2(:,:,els)*(intpd2_u.^3);
    else 
      cubic_term(:,els) = nek_mass(:,:,els)*(un(:,els).^3);  
    end
    b_cubic1 = extk(2)*cubic_term(:,els);
    b_cubic2 = extk(3)*cubic_lag1(:,els);
    b_cubic3 = extk(4)*cubic_lag2(:,els);
    b_cubic = deltat*(b_cubic1 + b_cubic2 + b_cubic3);

%   Quintic term
    if ifdealias  
      quintic_intp = nek_intpm1d3(:,:,els)*un(:,els);
      quintic_term(:,els) = nek_intgd3(:,:,els)*(quintic_intp.^5);
    else  
      quintic_term(:,els) = nek_mass(:,:,els)*(un(:,els).^5);
    end  
    b_quintic1 = extk(2)*quintic_term(:,els);
    b_quintic2 = extk(3)*quintic_lag1(:,els);
    b_quintic3 = extk(4)*quintic_lag2(:,els);
    b_quintic = deltat*(b_quintic1 + b_quintic2 + b_quintic3);

%   Forcing term (impulse)
    % spatial extent
    impulse_amp=1;
    nu = 0;
    sigma = 1.0;
    xofst_imp = 100;
    f_imp = impulse_amp*normpdf(xgll(:,els)-xofst_imp,nu,sigma);
    % temporal extent
    nu_t = 0;
    sigma_t = 0.01;
    tofst = 0.25;        % which part of the cycle
    tcycle = mod(time,Tosc)/Tosc;
%    timp = tcycle-tofst;
    timp = time;
    tnorm = normpdf(timp,nu_t,sigma_t)/normpdf(0,nu_t,sigma_t);
    forcing_term(:,els) = nek_mass(:,:,els)*f_imp*tnorm;

    b_forcing1 = extk(2)*forcing_term(:,els);
    b_forcing2 = extk(3)*forcing_lag1(:,els);
    b_forcing3 = extk(4)*forcing_lag2(:,els);
    b_forcing = deltat*(b_forcing1 + b_forcing2 + b_forcing3);

%   sum
    b(:,els) = -b1 + b2 + b_source + b_cubic - b_quintic + b_forcing;

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

  source_lag2 = source_lag1;
  source_lag1 = source_term; 

  cubic_lag2 = cubic_lag1;
  cubic_lag1 = cubic_term; 

  quintic_lag2 = quintic_lag1;
  quintic_lag1 = quintic_term;

  forcing_lag2 = forcing_lag1;
  forcing_lag1 = forcing_term;
%%

  periodic = 0; 
  b = DSSUM(b,nels,periodic);

%  b(1,1)    = 2e-14 + 1e-14*rand(1);    % random noise O(1e-6)
  b(1,1)    = 0;
%  b(1,1)    = 0.1*sin(time);
  big_b(1)  = b(1,1); 

  gl_mass(1,:) = 0;
  gl_mass(1,1) = 1;             % strong bc
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
  if mod(istep,iostep)==0 || istep==1 
    figure(h3)
    plot(ax1,xgll,un, '-k', 'MarkerSize', 16);
%    set(ax1,'Color', 'none');
    legend(ax1,['$T/T_{osc}=' num2str(time/Tosc) '$'], 'Location', 'Best', 'Interpreter', 'latex')

    mu_temp = mu_x0 + mu_time_0t + mu_diss*mu_xdiss + mu_c*xgll + mu_time_ct*xgll + mu_a*mu_xabs + mu_time_at*mu_xabs + mu_time_p*(xgll-pitch_x0);
    plot(ax2,xgll,mu_temp, '--r', 'MarkerSize', 16);
    set(ax2,'XAxisLocation', 'Top')
    set(ax2,'YAxisLocation', 'Right')
%    set(ax2,'XDir', 'reverse')  
    set(ax2,'Color', 'none')
    ylabel(ax1,'A')
    ylabel(ax2,'$\mu$', 'Interpreter', 'latex') 
%    ylim(ax2,[3 4.5])
    axes(ax2)
    set(ax2,'YGrid', 'on')
    set(ax2,'XGrid', 'on')
    pause(0.1)
  end

  if (mod(istep,isave)==0)
    u_save = [u_save un(:)];
    t_save = [t_save time];

    if (ifrenorm)
      [umax ind] = max(un(:));
      u_renorm = [u_renorm un(:)/umax];
    end


  end       

  %dbstop in me_periodic at 335
  iperiod=time/Tosc;
  if (floor(1/surfupd*(iperiod-nperiods))>=1)
    figure(10)
    surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
    view(2)
    colorbar

    if (ifrenorm)
      figure(11)
      surf(xgll(:),t_save/Tosc,transpose(u_renorm), 'EdgeColor', 'none')
      view(2)
      colorbar
    end

    pause(2)
    nperiods=nperiods+surfupd;
  end      
         
  time_it = toc(t_it);
end

figure(10)
surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
colorbar
view(2)
pause(2)

close all

%sfname=['GZ_' 'Uo' num2str(U0) '_muo' num2str(mu0) '_mut_o' num2str(mut_0) '_mu1_' num2str(mu1) '_mut_conv' num2str(mut_conv) '_mu_abs' num2str(mu_abs) '_mut_abs' num2str(mut_abs) '_xofst' num2str(xofst) '.mat'];
%save(sfname)



