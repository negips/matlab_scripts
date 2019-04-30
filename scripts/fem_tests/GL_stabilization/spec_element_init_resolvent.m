% Small test simulation.

clear
clc
close all

%addpath 'templates/'
ifplot = 0;

N=4;
lx1 = N+1;
npts=lx1;
nels = 4;
Nd=ceil(1.5*N);
Nd2=ceil((4*N+3)/2);
Nd3=ceil((6*N+3)/2);
ifboyd=0;

dof = N*nels+1;
nnodes = nels+1;

d_start = 0;
d_end = 20;
d_len = abs(d_end-d_start);

uniform = 1;

if uniform
  el_nodes = linspace(d_start,d_end,nnodes);
  el_size = diff(el_nodes);
else
  el_nodes = [-1 0.1 1];
  el_size = diff(el_nodes);
end

xgll = zeros(npts,nels);
gl_mass = zeros(dof,dof);
gl_gradm1 = zeros(dof,dof);
gl_stiff = zeros(dof,dof);
gl_intpm1d = zeros(dof,dof);
gl_intpm1d2 = zeros(dof,dof);
gl_lapl = zeros(dof,dof);
gl_cbc = zeros(dof,dof);
gl_lpbc = zeros(dof,dof);

mu_x = zeros(npts,nels);

io = sqrt(-1);    % iota

% \nu(dA/dx)
U     = 2.0;
cu    = 0.0;
nu = U + 2i*cu;

% \gamma(d^2A/dx^2)
cd = -1.;
gamma = 1.0 + 1i*cd;

% Source term
mu00   = 0.52;
mu0    = (mu00 - cu^2);
mu2    = -0.01;

% step for dissipative region
diss_xs=200;
diss_xe=300;
diss_xrise=20;
diss_xfall=20;

for i=1:nels

     close all;

     xst = el_nodes(i);
     xen = el_nodes(i+1);
     [MASS DXM1 DXM1D RXM1 gradm1 lpx FORC x_coeff Dx w1m1 xm1 JACM1 JACM1D xm1d xm1d2 GLL2Dealias Dealias2GLL GLL2Dealias2 Dealias2GLL2 GLL2Dealias3 Dealias2GLL3 gradm1d gradm1d2 gradm1d3 intgd intgd2 intgd3] = MESem1D2(N,Nd,Nd2,Nd3,xst,xen,ifboyd,ifplot);


     El(i).MASS   = MASS; 
     El(i).GRADM1 = gradm1;
     El(i).CONV   = MASS*gradm1;                % mass matrix included
     El(i).LAPL   = lpx;                        % mass matrix included
     El(i).INTPM1D  = GLL2Dealias;
     El(i).INTPM1D2 = GLL2Dealias2;  
     El(i).INTPM1D3 = GLL2Dealias3;  

     El(i).GRADM1D  = gradm1d;
     El(i).GRADM1D2 = gradm1d2;
     El(i).GRADM1D3 = gradm1d3;

     El(i).INTGD  = intgd;                      % mass matrix included
     El(i).INTGD2 = intgd2;                     % mass matrix included;
     El(i).INTGD3 = intgd3;                     % mass matrix included;

     El(i).xm1    = xm1; 
     El(i).xm1d   = xm1d; 
     El(i).xm1d2  = xm1d2;
     El(i).mu     = mu0 - mu1*xm1;
     El(i).mud    = mu0 - mu1*xm1d;

     nek_mass(:,:,i)     = MASS;
     nek_gradm1(:,:,i)   = gradm1;
     nek_conv(:,:,i)     = MASS*gradm1;
     nek_intpm1d(:,:,i)  = GLL2Dealias;
     nek_intpm1d2(:,:,i) = GLL2Dealias2;
     nek_intpm1d3(:,:,i) = GLL2Dealias3;
     nek_intgd(:,:,i)    = intgd;
     nek_intgd2(:,:,i)   = intgd2;
     nek_intgd3(:,:,i)   = intgd3;

     nek_lp(:,:,i)       = lpx;
     nek_mu(:,:,i)       = mu0*eye(length(xm1));
     nek_mud(:,:,i)      = mu0*eye(length(xm1d));

%    Time and space varying convective region      
     nek_mu_conv(:,:,i)  = diag(xm1);

%    On over-integration grid 
     nek_mud_conv(:,:,i) = diag(xm1d);

%    Add dissipative region at the end 
     s1                  =  smoothstep(xm1,diss_xs,diss_xs+diss_xrise); 
     s2                  = -smoothstep(xm1,diss_xe,diss_xe+diss_xfall); 
     diss_x              = s1+s2;     
     nek_mu(:,:,i)       = nek_mu(:,:,i) + mu_diss*diag(diss_x);
%    On over-integration grid 
     s1d                 =  smoothstep(xm1d,diss_xs,diss_xs+diss_xrise); 
     s2d                 = -smoothstep(xm1d,diss_xe,diss_xe+diss_xfall); 
     diss_xd             = s1d+s2d;     
     nek_mud(:,:,i)      = nek_mud(:,:,i) + mu_diss*diag(diss_xd);

     xgll(:,i) = xm1;

%    Just saving 
     s1                  =  smoothstep(xm1,step_xs,step_xs+step_xrise); 
     s2                  = -smoothstep(xm1,step_xe,step_xe+step_xfall);
     ds1                 =  smoothstep(xm1,diss_xs,diss_xs+diss_xrise); 
     ds2                 = -smoothstep(xm1,diss_xe,diss_xe+diss_xfall); 
     mu_x(:,i)           =  mu0*ones(length(xm1),1) + mu2/2*;
     mu_xdiss(:,i)       =  ds1 + ds2;
     mu_xabs(:,i)        =  s1 + s2; 

     gl_pos_j1 = (i-1)*N + 1;
     gl_pos_i1 = (i-1)*N + 1;

     gl_pos_j2 = (i)*N + 1;
     gl_pos_i2 = (i)*N + 1;
    
%     gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = MASS + gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
%     gl_gradm1(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = gradm1 + gl_gradm1(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
%     gl_intpm1d(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = GLL2Dealias + gl_intpm1d(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) ;
%     gl_intpm1d2(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = GLL2Dealias2 + gl_intpm1d2(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
%     gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = lpx + gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
%     gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = lpx + gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);

end

xall = unique(xgll);



