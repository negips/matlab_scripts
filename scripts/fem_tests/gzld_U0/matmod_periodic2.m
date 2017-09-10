% Small test simulation.

clear
%clc
close all

addpath 'templates/'

N=3;
npts = N+1;
nels = 3;
dof = N*nels+1;
nnodes = nels+1;
ifplot=0;

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
     el_nodes = [-1 0 1];
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

alt_eq0 = zeros(1,dof);
alt_eqn = zeros(1,dof);


nek_mass = zeros(npts,npts,nels);
nek_amass = zeros(npts,npts,nels);
nek_stiff = zeros(npts,npts,nels);
nek_conv = zeros(npts,npts,nels);
nek_conv_w = zeros(npts,npts,nels);
nek_forc = zeros(npts,npts,nels);
nek_lapl = zeros(npts,npts,nels);
nek_cbc = zeros(npts,npts,nels);
nek_lpbc = zeros(npts,npts,nels);

bord_mass      = zeros(npts,dof,2);
bord_amass     = zeros(npts,dof,2);
bord_stiff     = zeros(npts,dof,2);
bord_conv      = zeros(npts,dof,2);
bord_conv_w    = zeros(npts,dof,2);
bord_forc      = zeros(npts,dof,2);
bord_lapl      = zeros(npts,dof,2);
bord_cbc       = zeros(npts,dof,2);
bord_lpbc      = zeros(npts,dof,2);


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

     if periodic
          if i==1 || i==nels
               ind = floor(i/nels)+1;
               bord_mass      (:,(i-1)*N+1:i*N+1,ind)    = MASS;
               bord_amass     (:,(i-1)*N+1:i*N+1,ind)    = AMASS;
               bord_conv      (:,(i-1)*N+1:i*N+1,ind)    = CONV1;
               bord_conv_w    (:,(i-1)*N+1:i*N+1,ind)    = CONV_W;
               bord_forc      (:,(i-1)*N+1:i*N+1,ind)    = FORC;
               bord_lapl      (:,(i-1)*N+1:i*N+1,ind)    = LAPL;
               bord_cbc       (:,(i-1)*N+1:i*N+1,ind)    = CBC;
               bord_lpbc      (:,(i-1)*N+1:i*N+1,ind)    = LPBC;
          end
     end

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

xall=[];
for nel=1:nels
if nel==1
     ind=1:npts;
else
     ind=2:npts;
end
     xall=[xall;xgll(ind,nel)];
end

%%    Modify matrices for peridicity
if (periodic) && nels>2

     bord_conv(:,end,1)  = bord_conv(:,1,1);
     bord_conv(:,1,1)    = zeros(npts,1);
     bord_conv(:,:,1)    = flipud(bord_conv(:,:,1));
     bord_conv(:,1,2)    = bord_conv(:,end,2);
     bord_conv(:,end,2)  = zeros(npts,1);
     bord_conv(:,:,2)    = flipud(bord_conv(:,:,2));

     gl_conv(1:npts,:) = gl_conv(1:npts,:) + bord_conv(:,:,2);
     gl_conv(end-npts+1:end,:) = gl_conv(end-npts+1:end,:) + bord_conv(:,:,1);
     gl_conv=gl_conv(1:end-1,1:end-1);

%     bord_mass(:,end,1)  = bord_mass(:,1,1);
%     bord_mass(:,1,1)    = zeros(npts,1);
%     bord_mass(:,:,1)    = flipud(bord_mass(:,:,1));
%     bord_mass(:,1,2)    = bord_mass(:,end,2);
%     bord_mass(:,end,2)  = zeros(npts,1);
%     bord_mass(:,:,2)    = flipud(bord_mass(:,:,2));
%
%     gl_mass(1:npts,:) = gl_mass(1:npts,:) + bord_mass(:,:,2);
%     gl_mass(end-npts+1:end,:) = gl_mass(end-npts+1:end,:) + bord_mass(:,:,1);
     gl_mass=gl_mass(1:end-1,1:end-1);

%     gl_mass(:,1) = gl_mass(:,1) + gl_mass(:,end);
%     gl_mass(1,:) = gl_mass(1,:) + gl_mass(end,:);
%     gl_mass=gl_mass(1:end-1,1:end-1);
%
%     gl_amass(:,1) = gl_amass(:,1) + gl_amass(:,end);
%     gl_amass(1,:) = gl_amass(1,:) + gl_amass(end,:);
%     gl_amass=gl_amass(1:end-1,1:end-1);

     bord_forc(:,end,1)  = bord_forc(:,1,1);
     bord_forc(:,1,1)    = zeros(npts,1);
     bord_forc(:,:,1)    = flipud(bord_forc(:,:,1));
     bord_forc(:,1,2)    = bord_forc(:,end,2);
     bord_forc(:,end,2)  = zeros(npts,1);
     bord_forc(:,:,2)    = flipud(bord_forc(:,:,2));

     gl_forc(1:npts,:) = gl_forc(1:npts,:) + bord_forc(:,:,2);
     gl_forc(end-npts+1:end,:) = gl_forc(end-npts+1:end,:) + bord_forc(:,:,1);
     gl_forc=gl_forc(1:end-1,1:end-1);


%     gl_forc(:,1) = gl_forc(:,1) + gl_forc(:,end);
%     gl_forc(1,:) = gl_forc(1,:) + gl_forc(end,:);
%     gl_forc=gl_forc(1:end-1,1:end-1);
%
%     gl_lapl(:,1) = gl_lapl(:,1) + gl_lapl(:,end);
%     gl_lapl(1,:) = gl_lapl(1,:) + gl_lapl(end,:);
%     gl_lapl=gl_lapl(1:end-1,1:end-1);

     x=xall(1:end-1);
end
%% End of matrix modification


%clearvars -except N dt dt2 dc dc2 bc df_mat conv forc x abcde D

bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];

x2 = linspace(d_start,d_end,interp_pts);

% use a normalized gaussian
%mu = 0;
%sigma = 0.2;
%u0 = normpdf(x,mu,sigma)
%u0 = u0/max(u0);
%u0_int = normpdf(x2,mu,sigma);
%u0_int = u0_int/max(u0_int);

u0 = 0.1*sin(pi*x);
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

% Transformations due to periodicity
%if periodic

deltat = 0.005;
istep = 0;
nsteps = 10000;

ulag1 = 0*un;
ulag2 = 0*un;

ualag1 = 0*un;
ualag2 = 0*un;

%stiff = dt;
%astiff = dt2;

conv_s = 1.0;
chi = -0.0;
gl_dy_stiff = chi*gl_forc - conv_s*gl_conv;
%nek_dy_stiff = chi*nek_forc - conv_s*nek_conv_w;
%nek_bc = nek_cbc;

e1 = eig(inv(gl_mass)*gl_dy_stiff);
lambda_deltat = e1*deltat
clines = load('bdfk-neutral-curve.mat');

lambdar = real(e1)*deltat;
lambdai = imag(e1)*deltat;

h2 = figure;
plot(lambdar,lambdai, '*', 'MarkerSize', 16)
hold on
plot(clines.cline3(1,2:end),clines.cline3(2,2:end), 'r')
%xlim([min(lambdar) max(lambdar)]);
%ylim([min(lambdai) max(lambdai)]);

soln = 0*un;

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

     % Form RHS
     mass = gl_mass*bdfk(1);
     stiff = gl_dy_stiff;

     b11 = bdfk(2)*un; 
     b12 = bdfk(3)*ulag1; 
     b13 = bdfk(4)*ulag2;
     b1 = gl_mass*(b11+b12+b13);

     b21 = extk(2)*un;
     b22 = extk(3)*ulag1;
     b23 = extk(4)*ulag2;
     b2 = deltat*stiff*(b21 + b22 + b23);

     b = -b1 + b2;

     % set time integration operators

     restrt=N*nels;
     max_it = N*nels;
     tol = 1e-6;

     [soln cflag relres]  = gmres(mass,b,restrt,tol,max_it);

     display(['Relative Residual:' num2str(relres)])
     %soln = gl_mass\b;

     ulag2 = ulag1;
     ulag1 = un;
     un = soln;

     un2 = [un; un(1)];

     %% Just plot Points
     figure(h3)
     plot(xall,un2, 'ok', 'MarkerSize', 16);
     ylim([-0.15 0.15]);
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
     figure(h3)
     plot(x2,u_a,'b')
     ylim([-0.15 0.15]);
     hold off
     pause(0.002)
     %xaold = xanew;


end



