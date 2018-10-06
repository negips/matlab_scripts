clear all
close all
clc
% addpath ~/ProjectRE/Codes/matlab/wing_tools/
% addpath /scratch/negi/code_test/test_cases/wing_tools_ard/

path = './';
file   = 'saab_wing0.f00001';
filename = [path file];

Re_ref = 7.5e5;
ReC  = 7.5e05;       % chord-Re
scale_ratio = ReC/Re_ref;
C    = 1.;         % long-chord length [m]
nu   = 1.e-05;     % kin. viscosity
Q8   = ReC*nu/C;
rho8 = 1.;
pref = 4.4019e4;
nelz = 100;
lz   = 0.2;

casename = 'saab_wing';
%% read the mesh
[b,hdr,tag,N1,nel,Lvar,sts] = readnek_mpi(3,'le',filename,'xup');
N=N1-1;

% Rotation matrix for gll points
theta = 0.*pi/180;
rot = [cos(theta), -sin(theta); sin(theta) cos(theta)];

%Loop over all elements
%plotgll = zeros(length(b)*(N+1)^3,2);
for i = 1:length(b)
    vec_x = b{i}(:,:,:,1);
    vec_y = b{i}(:,:,:,2);
    vec_z = b{i}(:,:,:,3);
    
    vec_x_new = rot(1,1)*vec_x + rot(1,2)*vec_y;
    vec_y_new = rot(2,1)*vec_x + rot(2,2)*vec_y;
    vec_z_new = vec_z;
    
 %   EL(i).GLL_o(:,1) = squeeze(reshape(b{i}(:,:,:,1),N1^3,1,1)); % original data
 %   EL(i).GLL_o(:,2) = squeeze(reshape(b{i}(:,:,:,2),N1^3,1,1)); % original data
 %   EL(i).GLL_o(:,3) = squeeze(reshape(b{i}(:,:,:,3),N1^3,1,1)); % original data
    
    EL(i).GLL(:,1) = squeeze(reshape(vec_x_new,N1^3,1,1));
    EL(i).GLL(:,2) = squeeze(reshape(vec_y_new,N1^3,1,1));
    EL(i).GLL(:,3) = squeeze(reshape(vec_z_new,N1^3,1,1));
    
%    plotgll((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL(i).GLL;
end
%plot(plotgll(:,1),plotgll(:,2),'*','Markersize',2)
Nelm = length(b);
Ngp = N1;


%% Reading RANS Solution
disp('Write user-defined boundary conditions')
%load rans solution
%[xr,yr,ur,vr,wr,pr] = getrans(swa,'./solution_AoA_5.mat',lmult);
% solution = load('./sstkw-tr-aoa67.mat');
% solution = load('./sstkw-tr-aoa67_k05_a13.mat');
% solution = load('./sstkw-tr-aoa44_re400k.mat');
solution = load('./k-eps-aoa54_re750k.mat');

xr = solution.x;
yr = solution.y;

% Rotation matrix for rans grid and velocity
% positive anticlockwise
angle = -5.4;
disp(['Rotation angle: ' num2str(angle)])
theta = angle*pi/180;
rot = [cos(theta), -sin(theta); sin(theta) cos(theta)];

axis_x0 = 0.35;
axis_y0 = 0.034;
xrnew = xr-axis_x0;
yrnew = yr-axis_y0;

ind_min = find(xr==min(xr));

ur = solution.u;
vr = solution.v;
pr = solution.pressure;
wr = 0*ur;

l1 = length(xrnew);

for i = 1:l1

      coords = rot*[xrnew(i); yrnew(i)] + [axis_x0; axis_y0];
      xrnew(i) = coords(1);
      yrnew(i) = coords(2);

      vels = rot*[ur(i); vr(i)];
      ur(i) = vels(1);
      vr(i) = vels(2);
end

xr = xrnew;
yr = yrnew;

uinf = ur(ind_min);
pinf = pr(ind_min);

%normalize u
ur = ur/uinf;
vr = vr/uinf;
wr = wr/uinf;
%pr = pr/pinf;
%ric = solution(:,7);

% Normalize Ardeshir's pressure
facp8 = (rho8*Q8*Q8); % Wind-tunnel coords. -> plate coords.
pr = (pr-pinf)/facp8;

%Want positive W
wr = -wr;


%% write the initial condition
fpath = './';
filename_DNS = 'saab_wing.IC';
%writeic_2D_dns_sp
%writeic_2D_sp
writeic_sp
% writeic_sp_test

