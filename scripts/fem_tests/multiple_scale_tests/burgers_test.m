% Small test simulation.

clear
clc
close all

% Initialize 1D spectral element
spec_element_init_ms

A0 = 0.0;
A2 = 0.0;
A  = 1.0;
L  = d_len;

% Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];

x=xall;

%u0 = rand([length(xall),1]);
u0 = A0*ones([length(xall),1]);
% use a normalized gaussian
%nu = 0;
%sigma = 0.2;
%xofst = 5;
%u0 = 2*normpdf(xgll-xofst,nu,sigma);
%u0 = 1e-6*u0/max(max(u0));

periodic = 0; 
un = gs_scatter(u0,N,nels,periodic);


ulag1 = 0*un;
ulag2 = 0*un;

source_lag1 = 0*un;
source_lag2 = 0*un;

cubic_term = 0*un;
cubic_lag1 = 0*un;
cubic_lag2 = 0*un;

quintic_term = 0*un;
quintic_lag1 = 0*un;
quintic_lag2 = 0*un;

forcing_lag1 = 0*un;
forcing_lag2 = 0*un;

nlconv_term = 0*un;
nlconv_lag1 = 0*un;
nlconv_lag2 = 0*un;


nek_linear_stiff = -nek_conv - nek_lp;
%nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_linear_stiff;

ifcubic   = 0;    % Cubic term
ifquintic = 0;    % Quintic term
ifnlconv  = 0;    % Non-linear Convection

% Parameters
deltat = 0.00001;
istep = 0;
time = 0;
iostep = 1000;
isave  = 20;
OMEGA=0.50;
Tosc=2*pi/OMEGA;
nsteps=ceil(4*Tosc/deltat);
%nsteps=100000;
surfupd = 1000;               % make surface every surfupd of a period
verbose=1;                    % if verbose
verbosestep=5000;
ifpcg=0;
ifdealias=0;
ifrenorm=0;
ifsave = 0;
saveno = 0;

u_save   = [];
u_renorm = [];
t_save   = []; 

%% SOLVE
h3=figure;
ax1=axes;
axes(ax1);

nperiods=0;
for i = 1:nsteps
  t_it=tic;
  if (istep>0 && verbose && (mod(istep,verbosestep)==0))
    if (ifpcg)    
      disp(['Step, Time, Relative residual,iter,timing: ', num2str(istep), ', ' num2str(time), ', ' ...
        num2str(relres) ', ' num2str(pcgiter) ', ' num2str(time_it)]);
    else
      disp(['Step, Time, timing: ', num2str(istep,7), ', ' num2str(time,6), ', ' ...
        num2str(time_it,6)]);
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
    b_bdf = nek_mass(:,:,els)*(b11+b12+b13);

%   Conv + dissipation
    stiff(:,:,els) = nek_linear_stiff(:,:,els);
    b21 = extk(2)*un(:,els);
    b22 = extk(3)*ulag1(:,els);
    b23 = extk(4)*ulag2(:,els);
    b_stiff = deltat*stiff(:,:,els)*(b21 + b22 + b23);

%   Constant Source term
    if ifdealias  
      mu_term(:,els) = nek_intgd(:,:,els)*nek_mud(:,:,els)*nek_intpm1d(:,:,els)*un(:,els);
    else  
      mu_term(:,els) = nek_mass(:,:,els)*nek_mu(:,:,els)*un(:,els);
    end  

    source_term(:,els) = mu_term(:,els); 
    b_source = deltat*(extk(2)*source_term(:,els) + extk(3)*source_lag1(:,els) + extk(4)*source_lag2(:,els));


%   Cubic term
    if ifcubic
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
    else
      b_cubic = 0*un(:,els);
    end        

%   Quintic term
    if ifquintic
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
    else
      b_quintic = 0*un(:,els);
    end        

%   Nonlinear Convection
    
    if ifnlconv
      if 1 % ifdealias  
        nlconv_intp    = nek_intpm1d(:,:,els)*un(:,els);
        nlconv_gradm1d = nek_gradm1d(:,:,els)*un(:,els);

        nlconv_term(:,els) = nek_intgd(:,:,els)*(nlconv_intp.*nlconv_gradm1d);
      else
        nlconv_gradm1 = nek_gradm1(:,:,els)*un(:,els);
        nlconv_term(:,els) = nek_mass(:,:,els)*(un(:,els).*nlconv_gradm1);
      end  
      b_nlconv1 = extk(2)*nlconv_term(:,els);
      b_nlconv2 = extk(3)*nlconv_lag1(:,els);
      b_nlconv3 = extk(4)*nlconv_lag2(:,els);
      b_nlconv = deltat*(b_nlconv1 + b_nlconv2 + b_nlconv3);
    else
      b_nlconv = 0*un(:,els);
    end        


%   Forcing term (Initial impulse)
    if (istep==1)
      amp_imp= 0.;
      nu_imp = 0.;
      sigma_imp = 1.0;
      xofst_imp = -0;
      f_imp = amp_imp*normpdf(xgll(:,els)-xofst_imp,nu_imp,sigma_imp);
      forcing_term(:,els) = nek_mass(:,:,els)*f_imp;
    else
      forcing_term(:,els) = 0*un(:,els);
    end        

    b_forcing1 = extk(2)*forcing_term(:,els);
    b_forcing2 = extk(3)*forcing_lag1(:,els);
    b_forcing3 = extk(4)*forcing_lag2(:,els);
    b_forcing = -deltat*(b_forcing1 + b_forcing2 + b_forcing3);

%   sum for rhs
    b(:,els) = -b_bdf + b_stiff + b_source + b_cubic + b_quintic + b_nlconv + b_forcing;

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

  nlconv_lag2 = nlconv_lag1;
  nlconv_lag1 = nlconv_term;

  forcing_lag2 = forcing_lag1;
  forcing_lag1 = forcing_term;
%%

  periodic = 0; 
  b = DSSUM(b,nels,periodic);

% x=0 B.C.  
  b(1,1)    = A0 + A*sin(OMEGA*time);
  big_b(1)  = b(1,1); 

  gl_mass(1,:) = 0;
  gl_mass(1,1) = 1;             % strong bc

%% x=h B.C.
%  b(end,end) = A2;
%  big_b(end) = b(end,end); 
%
%  gl_mass(end,:)   = 0;
%  gl_mass(end,end) = 1;             % strong bc
 

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
    plot(ax1,xgll,un, '-b', 'MarkerSize', 16); hold on

%    ms_soln0 = ((A2-A0)/L - A/L*sin(OMEGA*time))*xgll + A0 + A*sin(OMEGA*time);
%    ms_soln1 = (A*L^2/nu)*cos(OMEGA*time)*(-1/3*(xgll/L) + 1/2*(xgll/L).^2 - 1/6*(xgll/L).^3);
%    ms_soln2 = -(A*L^4/nu^2)*sin(OMEGA*time)*(1/45*(xgll/L) - 1/18*(xgll/L).^3 + 1/24*(xgll/L).^4  - 1/120*(xgll/L).^5);
%    ms_soln3 = (A*L^6/nu^3)*cos(OMEGA*time)*(2/945*(xgll/L) - 1/270*(xgll/L).^3 + 1/360*(xgll/L).^5 - 1/720*(xgll/L).^6  + 1/5040*(xgll/L).^7);

    ms_soln0 = A*sin(OMEGA*time)*ones(size(xgll));      
    ms_soln1 = (A*L^2/nu)*cos(OMEGA*time)*(-xgll/L + 1/2*(xgll/L).^2);
    ms_soln2 = -(A*L^4/nu^3)*sin(OMEGA*time)*(1/3*xgll/L - 1/6*(xgll/L).^3 + 1/24*(xgll/L).^4);

    ms0      = ms_soln0;
    ms1      = ms0 + OMEGA*ms_soln1;
    ms2      = ms1 + (OMEGA^2)*ms_soln2;
  
%    plot(ax1,xgll,ms0, '--r', 'MarkerSize', 16); hold on
%    plot(ax1,xgll,ms1, '--g', 'MarkerSize', 16); hold on
    plot(ax1,xgll,ms2, '--m', 'MarkerSize', 16); hold on
%    plot(ax1,xgll,ms3, '--k', 'MarkerSize', 16); hold on
    ylim([-1. 1.])


%    legend(ax1,{'Numerical', '$0^{th}$','$1^{st}$','$2^{nd}$'}, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 12)
    legend(ax1,{'Numerical', 'Theoretical'}, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 12)
   
    ylabel(ax1,'$A$')
    xlabel(ax1,'$x$')
    pause(0.001)
    hold off

    if (ifsave)
      saveno = saveno+1; 
      figname = sprintf('multiple_scale_test%4.4d.png',saveno); 
      SaveFig(gcf,figname,'./plots',1)
    end        
    
  end

  if (mod(istep,isave)==0)
    u_save = [u_save un(:)];
    t_save = [t_save time];

    if (ifrenorm)
      [umax ind] = max(un(:));
      u_renorm = [u_renorm un(:)/umax];
    end


  end       
         
  time_it = toc(t_it);
end


%close all

%sfname=['GZ_' 'Uo' num2str(U0) '_muo' num2str(mu0) '_mut_o' num2str(mut_0) '_mu1_' num2str(mu1) '_mut_conv' num2str(mut_conv) '_mu_abs' num2str(mu_abs) '_mut_abs' num2str(mut_abs) '_xofst' num2str(xofst) '.mat'];
%save(sfname)



