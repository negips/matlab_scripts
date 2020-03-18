% Small test simulation.

clear
clc
close all

%addpath '../../'

tstart = tic;

display(['Time started: ' datestr(clock)])

SIZE           %    Polynomial order
re2            %    domain/decomposition/mapping/boundary conditions.

Re=1e+3;
nu=0; %1/Re;
ifboyd = 0;
destn  = 'plots/';
ifplot = 0;
ifdiv  = 1;
gltot  = 0;     %    Global number of elements 

for ii=1:nelv

     % clc 
     display(['Building Matrices for Element ', num2str(ii)])

     [mass nek_mass DXM1 DYM1 DXM1D DYM1D RXM1 RYM1 SXM1 SYM1 convx convy convall convxd convyd convalld convxd_new convyd_new convalld_new Cfx Cfy gradm1x gradm1y gradm1xd gradm1yd intpm1d intpd2m1 wtsvecd massd nek_conv lpx lpy lpall nek_lp lpbc forc NxNy_nodal2spec NxdNyd_spec2nodal Nx_GLL2Gauss x_coeff y_coeff Dx Dy w2m1 xm1 ym1 JACM1 JACM1D xm1d ym1d] = MESem2D8(Nx,Ny,Nxd,Nyd,El(ii).xc,El(ii).yc,ifboyd,ifplot);

     El(ii).mass = mass;      
     El(ii).nek_mass = nek_mass;
     El(ii).Dx = Dx;
     El(ii).Dy = Dy;
     El(ii).x_coeff = x_coeff;
     El(ii).y_coeff = y_coeff;
     El(ii).DXM1 = DXM1;
     El(ii).DYM1 = DYM1;
     El(ii).DXM1D = DXM1D;
     El(ii).DYM1D = DYM1D;
     El(ii).RXM1 = RXM1;
     El(ii).RYM1 = RYM1;
     El(ii).SXM1 = SXM1;
     El(ii).SYM1 = SYM1;
     El(ii).convx = convx;
     El(ii).convy = convy;
     El(ii).convall = convall;
     El(ii).forc = forc;
     El(ii).NxNy_nodal2spec = NxNy_nodal2spec;
     El(ii).NxdNyd_spec2nodal = NxdNyd_spec2nodal;
     El(ii).Nx_GLL2Gauss = Nx_GLL2Gauss;

     El(ii).convxd = convxd;
     El(ii).convyd = convyd;
     El(ii).convalld = convalld;

     El(ii).convxd_new = convxd_new;
     El(ii).convyd_new = convyd_new;
     El(ii).convalld_new = convalld_new;

     El(ii).Cfx = Cfx;
     El(ii).Cfy = Cfy;

     El(ii).gradm1x = gradm1x;
     El(ii).gradm1y = gradm1y;
     El(ii).gradm1xd = gradm1xd;
     El(ii).gradm1yd = gradm1yd;
     El(ii).intpm1d = intpm1d;
     El(ii).intpd2m1 = intpd2m1;
     El(ii).wtsvecd = wtsvecd;
     El(ii).massd = massd;
     El(ii).xm1d = xm1d;
     El(ii).ym1d = ym1d;
     El(ii).JACM1D = JACM1D;

     El(ii).nek_conv = nek_conv;
     El(ii).laplx = lpx;
     El(ii).laply = lpy;
     El(ii).lapl = lpall;
     El(ii).nek_lp = nek_lp;
     El(ii).xm1 = xm1;
     El(ii).ym1 = ym1;
     El(ii).JACM1 = JACM1;
     El(ii).w2m1 = w2m1;

     El(ii).un = 0*xm1;
     El(ii).ulag1 = El(ii).un;
     El(ii).ulag2 = El(ii).un;
     El(ii).soln = El(ii).un;
     El(ii).b = El(ii).un;
     El(ii).bl = El(ii).un;
     El(ii).scrtch1 = El(ii).un;
     El(ii).scrtch2 = El(ii).un;

     %% Calculate divergence
     div = El(ii).gradm1x*Cfx(:) + El(ii).gradm1y*Cfy(:);
     El(ii).div = reshape(div,Nx+1,Ny+1);
     divd = El(ii).gradm1xd*Cfx(:) + El(ii).gradm1yd*Cfy(:);
     El(ii).divd = reshape(divd,Nxd+1,Nyd+1);

%    Only for testing. Delete 
%     convection = El(ii).convalld*El(ii).Cfx(:); 
%     El(ii).convection = reshape(convection,Nx+1,Ny+1);

     %%
end

return

%    Build multiplicity
for ii = 1:nelv
     [El(ii).lglmap gltot] = GlobalNo(ii,Nx,Ny,nelx,nely,xperiodic,yperiodic,gltot,El);
end
El = vmult(El,nelv,lx1,ly1);

%    inidices/Global nos of boundary nodes
El = BdryIndex(El,nelv,lx1,ly1);

% initial Condition
El = SetICs(El,nelv);

% Check Domain volume
volume = GetVolume(El,nelv);
disp(['Domain Volume: ' num2str(volume)])

display(['Initialization completed: ' datestr(clock)])
toc(tstart)

%------------------------------------------------------------ 
% Time Stepper

% Initialize
x0 = 0*El.un(:);
v0 = 0*x0;
a0 = 0*x0;















