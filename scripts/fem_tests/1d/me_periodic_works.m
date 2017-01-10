% Small test simulation.

clear
clc
close all

addpath 'templates/'
ifplot = 0;

N=2;
npts = N+1;
nels = 2;

dof = N*nels+1;
nnodes = nels+1;

d_start = -1;
d_end = 1;
d_len = abs(d_end-d_start);

periodic = 1;

interp_pts = 200;

uniform = 1;

if uniform
     el_nodes = linspace(d_start,d_end,nnodes);
     el_size = diff(el_nodes);
else
%     el_nodes = lglnodes(nels);
%     el_nodes = el_nodes(end:-1:1);
     el_nodes = [-1 0.1 1];
     el_size = diff(el_nodes);
end

xgll = zeros(npts,nels);
gl_mass = zeros(dof,dof);
gl_amass = zeros(dof,dof);
gl_stiff = zeros(dof,dof);
gl_conv = zeros(dof,dof);
gl_forc = zeros(dof,dof);
gl_lapl = zeros(dof,dof);
gl_cbc = zeros(dof,dof);
gl_lpbc = zeros(dof,dof);

nek_mass = zeros(npts,npts,nels);
nek_amass = zeros(npts,npts,nels);
nek_stiff = zeros(npts,npts,nels);
nek_conv = zeros(npts,npts,nels);
nek_conv_w = zeros(npts,npts,nels);
nek_forc = zeros(npts,npts,nels);
nek_lapl = zeros(npts,npts,nels);
nek_cbc = zeros(npts,npts,nels);
nek_lpbc = zeros(npts,npts,nels);

for i=1:nels

     close all;

     xst = el_nodes(i);
     xen = el_nodes(i+1);
     [MASS AMASS CONV_W CONV1 LAPL CBC LPBC FORC x_p_coeff dx_p_coeff xorg wts] = MEFem(N,xst,xen,ifplot);
     xgll(:,i) = el_nodes(i) + (xorg+1)/2*el_size(i);

     nek_mass(:,:,i) = MASS;
     nek_amass(:,:,i) = AMASS;
     nek_conv(:,:,i) = CONV1;
     nek_conv_w(:,:,i) = CONV_W;
     nek_forc(:,:,i) = FORC;
     nek_lapl(:,:,i) = LAPL;
     nek_cbc(:,:,i) = CBC;
     nek_lpbc(:,:,i) = LPBC;

     gl_pos_j1 = (i-1)*N + 1;
     gl_pos_i1 = (i-1)*N + 1;

     gl_pos_j2 = (i)*N + 1;
     gl_pos_i2 = (i)*N + 1;
    
     gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = MASS + gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
     gl_amass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = AMASS + gl_amass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
     gl_conv(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = CONV1 + gl_conv(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) ;
     gl_forc(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = FORC + gl_forc(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
     gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = LAPL + gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);

end
break
xall = unique(xgll(:));

%clearvars -except N dt dt2 dc dc2 bc df_mat conv forc x abcde D

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
%mu = 0;
%sigma = 0.2;
%u0 = normpdf(x,mu,sigma)
%u0 = u0/max(u0);
%u0_int = normpdf(x2,mu,sigma);
%u0_int = u0_int/max(u0_int);

u0 = 0.1*sin(pi*xgll);
u0_int = 0.1*sin(pi*x2);

%u0 = 1 - x.^2;
%u0_int = 1 - x2.^2;

un = u0;

%% testing
%xtemp = transpose(linspace(-1,1,200));
%l2 = length(x2);
%basis = zeros(l2,N+1);
%for i = 0:N 
%     val = zeros(l2,1);
%for j = 0:N
%     val = val + un(i+1)*x_p_coeff(j+1,i+1)*x2.^j;
%end
%     basis(:,i+1) = val;
%end
%
%vals = sum(basis,2);
%h1=figure;
%plot(x2,vals);
%hold on
h_isoln = figure;
plot(x2,u0_int,'o');
title('initial solution');

%% testing plot derivative
%xtemp = transpose(linspace(-1,1,200));
%l2 = length(x2);
%basis = zeros(l2,N+1);
%for i = 0:N
%val = zeros(l2,1);
%for j = 1:N              % first index val is zero
%     val = val + un(i+1)*dx_p_coeff(j+1,i+1)*x2.^(j-1);
%end
%basis(:,i+1) = val;
%end
%vals = sum(basis,2);

u0_deriv = pi*0.1*cos(pi*x2);
%u0_deriv = -2*x2;

%h_deriv=figure;
%plot(x2,u0_deriv,'ok');
%title('initial derivatives')

deltat = 0.002;
istep = 0;
nsteps = 10000;

ulag1 = 0*un;
ulag2 = 0*un;

ualag1 = 0*un;
ualag2 = 0*un;

%stiff = dt;
%astiff = dt2;

conv_s = 2.0;
chi = -0.0;
re = 1e+1;
nek_dy_stiff = chi*nek_forc - conv_s*nek_conv;
nek_bc = nek_cbc;
mass = 0*nek_mass;
stiff = 0*nek_stiff;


%% Just to get eigenvalues:
     gl_conv(1,1) = 0;
     gl_conv(end,end) = 0;
     gl_dy_mass = gl_mass - 1/re*gl_lapl;
     gl_dy_stiff = chi*gl_forc - conv_s*gl_conv;
     e1 = eig(inv(gl_dy_mass)*gl_dy_stiff);
     lambda_deltat = e1*deltat
     clines = load('bdfk-neutral-curve.mat');

     lambdar = real(e1)*deltat;
     lambdai = imag(e1)*deltat;

     h2 = figure;
     plot(lambdar,lambdai,'*', 'MarkerSize', 16)
     hold on
     plot(clines.cline3(1,2:end),clines.cline3(2,2:end), 'r')
     %xlim([min(lambdar) max(lambdar)]);
     %ylim([min(lambdai) max(lambdai)]);
%%
break

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

          stiff(:,:,els) = nek_dy_stiff(:,:,els); 
          b21 = extk(2)*un(:,els);
          b22 = extk(3)*ulag1(:,els);
          b23 = extk(4)*ulag2(:,els);
          b2 = deltat*stiff(:,:,els)*(b21 + b22 + b23);

          b(:,els) = -b1 + b2;

     end

     b = DSSUM(b,nels,1);

     M=[];                         % preconditioner
     restrt=N*nels;
     max_it = 5;
     tol = 1e-9;

     [soln err iter] = me_solve(mass,un,b,max_it,tol);

     ulag2 = ulag1;
     ulag1 = un;
     un = soln;

     %% Just plot Points
     figure(h3)
     plot(xgll,un, 'ok', 'MarkerSize', 16);
     hold on
     %plot(x,una2, '*r', 'MarkerSize', 16);

     xanew = x2 - conv_s*time;
     for j=1:length(xanew)
          if xanew(j)>1
               while(xanew(j)>1)
                    xanew(j) = xanew(j)-2;
               end
          elseif xanew(j)<-1
               while (xanew(j)<-1)
                    xanew(j) = xanew(j) + 2;
               end
          end
     end

     %u_a = 1-xanew.^2;
     u_a = 0.1*sin(pi*xanew);

     plot(x2,u_a,'b')
     ylim([-0.15 0.15]);
     hold off
     pause(0.005)
     %xaold = xanew;

     %dbstop in me_periodic at 335


end



