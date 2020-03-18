% Small test simulation.

clear
clc
close all

% Initialize 1D spectral element
spec_element_init_time

% Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];

%who

x=xall;

u0 = rand([length(xall),1]);
% use a normalized gaussian
k0 = 1.5;
nu = 0;
sigma = 2.5;
xofst = 50;
u0 = 2*normpdf(xgll-xofst,nu,sigma).*cos(k0*xgll);
u0 = 1e-0*u0/max(max(u0));
xf = INTP*xgll;
u0f= INTP*u0;
uen= envelope(u0f(:));

fig0 = figure(1);
plot(xf,u0f); hold on
plot(xf(:),uen, '--k')

%return

un = u0;

ulag1 = 0*un;
ulag2 = 0*un;

stiff_lag1  = 0*un;
stiff_lag2  = 0*un;

source_lag1 = 0*un;
source_lag2 = 0*un;

cubic_lag1 = 0*un;
cubic_lag2 = 0*un;

quintic_lag1 = 0*un;
quintic_lag2 = 0*un;

forcing_lag1 = 0*un;
forcing_lag2 = 0*un;

%conv_s = U0;
%nek_linear_stiff = -conv_s*nek_conv - nek_lp;
%nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_mass;

% Parameters
deltat = 0.0002;
istep = 0;
nsteps = 100000;
time = 0;
iostep = 500;
isave  = 100;
OMEGA=0.1;
Tosc=2*pi/OMEGA;
%nsteps=ceil(1.5*Tosc/deltat);
surfupd = 1000;               % make surface evry surfupd of a period
verbose=1;
verbosestep=500;
ifpcg=0;
ifdealias=0;
ifrenorm=0;

u_save = [];
u_renorm = [];
t_save = [];
pk_x   = [];
pk_u   = [];

%% SOLVE
h3=figure(2);
ax1=axes;

scale  = 1.0;
scfac  = 1.0e+5;
scfac2 = 1.0e-5;

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
    U0 = 2.0;
    dU = 0.0;
    gamma = 1.0;
    nek_linear_stiff = -(U0 + dU*sin(OMEGA*time))*nek_conv - gamma*nek_lp;

    stiff(:,:,els) = nek_linear_stiff(:,:,els);
    stiff_term(:,:,els) = stiff(:,:,els)*un(:,els);
    b21 = extk(2)*stiff_term(:,:,els);
    b22 = extk(3)*stiff_lag1(:,els);
    b23 = extk(4)*stiff_lag2(:,els);
    b2 = deltat*(b21 + b22 + b23);

%   Spatially constant Source term
    mu0 = 0.50;
    dmu = 0.25;
    mu  = mu0 + dmu*sin(OMEGA*time);
    if ifdealias  
      mu_term(:,els) = mu*nek_intgd(:,:,els)*nek_intpm1d(:,:,els)*un(:,els);
    else  
      mu_term(:,els) = mu*nek_mass(:,:,els)*un(:,els);
    end  

%   Space varying Source term    
    mux = -0.000;
    if ifdealias  
      mux_term(:,els) = mux*nek_intgd(:,:,els)*nek_muxd(:,:,els)*nek_intpm1d(:,:,els)*un(:,els);
    else  
      mux_term(:,els) = mux*nek_mass(:,:,els)*nek_mux(:,:,els)*un(:,els);
    end  

    mux2 = -0.000;
    if ifdealias  
      mux2_term(:,els) = mux2*nek_intgd(:,:,els)*nek_muxd2(:,:,els)*nek_intpm1d(:,:,els)*un(:,els);
    else  
      mux2_term(:,els) = mux2*nek_mass(:,:,els)*nek_mux2(:,:,els)*un(:,els);
    end

    source_term(:,els) = mu_term(:,els) + mux2_term(:,els) + mux2_term(:,els);

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
    b_cubic = 0*deltat*(b_cubic1 + b_cubic2 + b_cubic3);

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
    b_quintic = 0*deltat*(b_quintic1 + b_quintic2 + b_quintic3);

%   sum
    b(:,els) = -b1 + b2 + b_source + b_cubic - b_quintic;

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

  stiff_lag2 = stiff_lag1;
  stiff_lag1 = stiff_term; 

  source_lag2 = source_lag1;
  source_lag1 = source_term; 

  cubic_lag2 = cubic_lag1;
  cubic_lag1 = cubic_term; 

  quintic_lag2 = quintic_lag1;
  quintic_lag1 = quintic_term;

%  forcing_lag2 = forcing_lag1;
%  forcing_lag1 = forcing_term;


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

  if (mod(istep,isave)==0)
    u_save = [u_save scale*un(:)];
    t_save = [t_save time];

    uint = INTP*un;
    fl = 50000;
    uen  = envelope(uint(:),fl,'analytic');

    [pks ind] = max(uen);
    pk_u = [pk_u scale*pks];

%    xint = INTP*xgll;
    pk_x = [pk_x xf(ind)];

%   Just plot Points
    if (mod(istep,iostep)==0)
      figure(h3)
      pl=plot(xf,uint, '-b', 'MarkerSize', 16); hold on
      pl2=plot(xf(:),uen, '--k', 'MarkerSize', 16); hold off
%      set(ax1,'Color', 'none');
      legend(pl(1),{['$T/T_{osc}=' num2str(time/Tosc) '$']}, 'Location', 'NorthEast', 'Interpreter', 'latex')

      ylabel('$A$')
      pause(0.01)
    end  
  end



  %dbstop in me_periodic at 335
  if (mod(istep,surfupd)==0 && istep>0)
    figure(10)
    surf(xgll(:),t_save/Tosc,log10(abs(transpose(u_save))+1.0e-15), 'EdgeColor', 'none')
    view(2)
    colorbar

    figure(11)
    plot(pk_x,t_save);
    grid on
    xlabel('$x$')
    ylabel('$t$')

    figure(12)
    semilogy(t_save,abs(pk_u)+1.0e-15);
    grid on
    xlabel('$t$')
    ylabel('$u^{pk}$')

    pause(0.01)
  end      
         
  time_it = toc(t_it);

  umax = max(abs(v));
  if (umax>scfac)
    scale = scale*scfac;

    ulag2 = ulag2/scfac;
    ulag1 = ulag1/scfac;

    stiff_lag2 = stiff_lag2/scfac;
    stiff_lag1 = stiff_lag1/scfac; 

    source_lag2 = source_lag2/scfac;
    source_lag1 = source_lag1/scfac; 

    cubic_lag2 = cubic_lag2/scfac;
    cubic_lag1 = cubic_lag1/scfac; 

    quintic_lag2 = quintic_lag2/scfac;
    quintic_lag1 = quintic_lag1/scfac;

    un = un/scfac;
  end        

  if (umax<scfac2)

    scale = scale*scfac2;

    ulag2 = ulag2/scfac2;
    ulag1 = ulag1/scfac2;

    stiff_lag2 = stiff_lag2/scfac2;
    stiff_lag1 = stiff_lag1/scfac2; 

    source_lag2 = source_lag2/scfac2;
    source_lag1 = source_lag1/scfac2; 

    cubic_lag2 = cubic_lag2/scfac2;
    cubic_lag1 = cubic_lag1/scfac2; 

    quintic_lag2 = quintic_lag2/scfac2;
    quintic_lag1 = quintic_lag1/scfac2;

    un = un/scfac2;
  end        


end

%figure(10)
%surf(xgll(:),t_save/Tosc,transpose(u_save), 'EdgeColor', 'none')
%colorbar
%view(2)
%pause(2)

%close all

dt = deltat*isave;
pk_v  = diff(pk_x)/dt;
pk_v2 = pk_v;
tsave2 = t_save(2:end);
nsave = length(pk_v);
tmp = 0*pk_v;

for i=1:nsave-20
  ps   = sort(pk_v(i:i+20));
  pkmi = mean(ps(5:15));
%  stdx = std(pk_v(i-9:i));
  stdx = std(ps(5:15));
  tmp(i) = stdx;
  tmp2(i)= pkmi;
  if pk_v(i)<(pkmi+5*stdx) && pk_v(i)>(pkmi-5*stdx)
    pk_v2(i)= pk_v(i);
  else
    pk_v2(i)=-100;
  end
end
ii = find(pk_v2~=-100);
Ut = (U0 + dU*sin(OMEGA*t_save));

figure
plot(tsave2(ii),pk_v2(ii), '--b'); hold on
plot(t_save,Ut,'k')

sfname=['GZ_time_k01.mat'];
save(sfname)


