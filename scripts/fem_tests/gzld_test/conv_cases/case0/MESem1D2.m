function [mass DXM1 DXM1D RXM1 gradm1 lpx forc x_coeff Dx w1m1 xm1 JACM1 JACM1D xm1d xm1d2 GLL2Dealias Dealias2GLL GLL2Dealias2 Dealias2GLL2 GLL2Dealias3 Dealias2GLL3 gradm1d gradm1d2 gradm1d3 intgd intgd2 intgd3] = MESem1D2(Nx,Nxd,Nxd2,Nxd3,xs,xe,ifboyd,ifplot)

%addpath '../../';

%    Definition of basis functions:
%    Definition of nodal values

%    U0-----U1-----U2

% For full matrix:
% Vector: {U0,
%          U1,
%          U2}

tic

epstol=1e-30;

if Nx>0
     [x wx p]= lglnodes(Nx);
     x =x(end:-1:1);
     wx =wx(end:-1:1);
else
     x=0;
     wx=1;
end

w1m1 = zeros(Nx+1,1);
for i=0:Nx
  w1m1(i+1) = wx(i+1);
end

N = Nx;                  % temporary

% Build polynomial coefficients for x 
A1 = [];
for i = 0:Nx
     A1= [A1 x.^i];
end
x_coeff = zeros(Nx+1);
%
%for i = 0:Nx
%  b = zeros(Nx+1,1);
%  b(i+1) = 1;
%
%  solns = A1\b;
%  x_coeff(:,i+1) = solns;
%end
%---------------------------------------- 
% Build derivative matrix (d/dx)
Dx = [];
%for i = 0:Nx
%  j = (i-1);
%  t1 = i*x_coeff(i+1,:);
%  nans = isnan(t1);
%  if max(nans)
%    ind = find(nans);
%    t1(ind) = 0;
%  end
%  Dx= [Dx; t1];
%end
%---------------------------------------- 
disp('Calculating Jacobian')
%% Derivative operator d/dx
DXM1 = zeros(Nx+1,Nx+1);
d0 = Nx*(Nx+1)/4;
for ii=0:Nx
  for jj=0:Nx
    DXM1(ii+1,jj+1) = 0;
    if ii~=jj
      lni = legendrePoly(Nx,x(ii+1));
      lnj = legendrePoly(Nx,x(jj+1));
      DXM1(ii+1,jj+1) = lni(end)/lnj(end)/(x(ii+1)-x(jj+1));
    end
    if ii==jj && ii==0
      DXM1(ii+1,jj+1) = -d0;
    end
    if ii==jj && ii==Nx
      DXM1(ii+1,jj+1) = d0;
    end
  end
end
%---------------------------------------- 
xm1 = getlgll1D(Nx,xs,xe);          % linear mapping

XRM1 = DXM1*xm1;
JACM1 = 2/(xe-xs);

%---------------------------------------- 
% All these factors need to be divided by the jacobian (point wise multiplication).
% It get multiplied by the jacobian during integral.
% Hence it is skipped right now.
RXM1 = 1./XRM1;

%---------------------------------------- 

% Mass matrix (in front of d/dt)
% This is a big but sparse matrix
disp('Calculating Mass matrix')
mass = zeros((Nx+1),(Nx+1));

for n=0:Nx
  mass(n+1,n+1) = JACM1*w1m1(n+1);
end    
%-------------------------------------------------- 

%% Over-integration matrices
if Nxd>1
  [xd wxd p]= lglnodes(Nxd);
  xd =xd(end:-1:1);
  wxd =wxd(end:-1:1);
else
  xd=[0];
  wxd=[1];
end

%---------------------------------------- 
%% Derivative operator d/dx
DXM1D = zeros(Nxd+1,Nx+1);
%----------------------------------------
xm1d = getlgll1D(Nxd,xs,xe);

JACM1D = 2/(xe-xs);
%-------------------------------------------------- 
%% Over-integration matrix 2
if Nxd2>1
  [xd2 wxd2 p]= lglnodes(Nxd2);
  xd2 =xd2(end:-1:1);
  wxd2 =wxd2(end:-1:1);
else
  xd2=[0];
  wxd2=[1];
end
%---------------------------------------- 
xm1d2 = getlgll1D(Nxd2,xs,xe);

JACM1D2 = 2/(xe-xs);
%---------------------------------------- 
% All these factors need to be divided by the jacobian (point wise multiplication).
% It get multiplied by the jacobian during integral.
% Hence it is skipped right now.

%-------------------------------------------------- 
%% Over-integration matrix 3
if Nxd3>1
  [xd3 wxd3 p]= lglnodes(Nxd3);
  xd3 =xd3(end:-1:1);
  wxd3 =wxd3(end:-1:1);
else
  xd3=[0];
  wxd3=[1];
end
%---------------------------------------- 
xm1d3 = getlgll1D(Nxd3,xs,xe);

JACM1D3 = 2/(xe-xs);
%---------------------------------------- 
% All these factors need to be divided by the jacobian (point wise multiplication).
% It get multiplied by the jacobian during integral.
% Hence it is skipped right now.

%% Dealiasing
%-------------------------------------------------- 
%disp('Calculating dealiased linear convective matrix')
% Interpolate Geometrical factors to over-integration grid
%-------------------------------------------------- 
I_JACM1D = JACM1D;


% Mass matrix for dealiased operator
%-------------------------------------------------- 
massd = zeros((Nxd+1),(Nxd+1));
for m=0:Nxd
  massd(m+1,m+1) = I_JACM1D*wxd(m+1);
end 

%-------------------------------------------------- 
massd2 = zeros((Nxd2+1),(Nxd2+1));
for m=0:Nxd2
  massd2(m+1,m+1) = I_JACM1D*wxd2(m+1);
end 
%-------------------------------------------------- 
massd3 = zeros((Nxd3+1),(Nxd3+1));
for m=0:Nxd3
  massd3(m+1,m+1) = I_JACM1D*wxd3(m+1);
end 
%-------------------------------------------------- 

% Just derivative operator.
% No integration
% du/dx (on the dealiased grid).

disp('Calculating derivative operators')
%-------------------------------------------------- 
gradm1 = zeros((Nx+1),(Nx+1));

%% GLL points
%for k =0:Nx                        % points along r
%  wtr = wx(k+1);
%  rk = x(k+1);
%
%  for i = 0:Nx
%    Li = x_coeff(:,i+1);          % Trial Function (ur)
%    dLi = Dx(:,i+1);
%
%%   Test function/derivative values.
%    ifderiv =0;
%    Li_rk = FuncEval(Li,rk,ifderiv);
%
%    ifderiv =1;
%    dLi_rk = FuncEval(dLi,rk,ifderiv);
%
%    deriv_val = dLi_rk*RXM1(k+1);
%
%    posx = k + 1;
%    posy = i + 1;
%%   Since we build this very generally...
%%   Numerical errors occur (O(1e-14))
%    if (abs(deriv_val)>epstol)  
%      gradm1(posx,posy) = deriv_val;
%    end  
%  end
%end

%% testing
%% Build legendre transform
pht = zeros(Nx+1,Nx+1);
dpht = zeros(Nx+1,Nx+1);
nx = Nx+1;
lx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=x(j);
  Lj = legendrePoly(n,z);
  dLj=legendreDeriv(n,z)*RXM1(j);
  pht(j,:)  = transpose(Lj);
  dpht(j,:) = transpose(dLj);
end

InvLegendreTransform = pht;             % Spectral to nodal 
LegendreTransform = inv(pht);           % Nodal to Spectral
gradm1 = dpht*LegendreTransform;    % gradm1

%dbstop in MESem1D at 331
%%

%% Laplacian term                  % Integration by parts:  dudx.dvdx
%---------------------------------------------------------------------- 
disp('Calculating laplacian')
lpx = zeros((Nx+1),(Nx+1));
trial_fcn = zeros((Nx+1),1);
deriv_trial = zeros((Nx+1),1);
trial_mat = zeros((Nx+1),(Nx+1));
geom_fac = diag(JACM1*w1m1);
for l=1:(Nx+1)         % trial function
  trial = trial_fcn;
  trial(l) = 1;
  deriv_trial = gradm1*trial;
  trial_mat(l,:) = transpose(deriv_trial);
end
lpx = trial_mat*geom_fac*gradm1;

%--------------------------------------------------
lpbc = zeros(N+1);             % Boundary terms for laplacian operator. vdudx(1) - vdudx(-1)




%% Build forcing
%---------------------------------------------------------------------- 
%% Build basis change matrix
%% Taken directly from nek
%% In the x-direction
nx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:nx
  z=x(j);
  Lj = legendrePoly(n,z);
  kj = kj+1;
  pht(kj) = Lj(1);
  kj = kj+1;
  pht(kj) = Lj(2);
  for k = 3:nx
    kj = kj+1;
    pht(kj) = Lj(k);% - Lj(k-2);
  end
end

pht = reshape(pht,nx,nx);
Nx_spec2nodal = transpose(pht);
Nx_nodal2spec = inv(Nx_spec2nodal);

% Filter function. Modes and amplitudes
alpha_x=1.0;
Gx = eye(Nx+1);
kc=Nx-2;
for k=1:Nx+1
  if (k>kc)    
    alpha_x=((k-kc)/(Nx+1-kc))^2;
  else
    alpha_x=0;
  end
  Gx(k,k)=1. - alpha_x;
end  

% Build matrix which applies filter on Boyd transformed basis.
% And then transforms back to physical space.
phi = pht;

if ifboyd
  for i = 3:nx
    for j = 1:nx
      phi(i,j) = pht(i,j) - pht(i-2,j);
    end
  end
end
boydstonodal_x = transpose(phi);
nodaltoboyds_x = inv(boydstonodal_x);

fil_mat_x = boydstonodal_x*Gx*nodaltoboyds_x;

forc = mass*(eye((Nx+1))-fil_mat_x);
%-------------------------------------------------- 

%% Build Over-integration matrices
% Nxd1
pht = zeros(Nx+1,Nxd+1);
nx = Nxd+1;
lx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=x(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end

dealiased_intpx = pht;              % Nxd spectral to Nx nodal
%% Build spectral transform matrix on Nxd
pht = zeros(Nxd+1,Nxd+1);
nx = Nxd+1;
lx = Nxd+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
Nxd_spec2nodal = pht;              % Nxd spectral to Nxd nodal
Nxd_nodal2spec = inv(pht);         % Nxd nodal to Nxd spectral

Dealias2GLL = dealiased_intpx*Nxd_nodal2spec;
intgd = Dealias2GLL*massd;
%% Spectral Nx to nodal Nxd
pht = zeros(Nxd+1,Nx+1);
nx = Nx+1;
lx = Nxd+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end

Nxspec_to_Nxdnodal = pht;              % Nx spectral to Nxd nodal

GLL2Dealias = Nxspec_to_Nxdnodal*LegendreTransform;
gradm1d = GLL2Dealias*gradm1;

%-------------------------------------------------- 
% Nxd2
pht = zeros(Nx+1,Nxd2+1);
nx = Nxd2+1;
lx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=x(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
dealiased_intpx2 = pht;            % Nxd2 spectral to Nx nodal
%% Build spectral transform matrix on Nxd2
pht = zeros(Nxd2+1,Nxd2+1);
nx = Nxd2+1;
lx = Nxd2+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd2(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
Nxd2_spec2nodal = pht;              % Nxd spectral to Nxd nodal
Nxd2_nodal2spec = inv(pht);         % Nxd nodal to Nxd spectral

Dealias2GLL2 = dealiased_intpx2*Nxd2_nodal2spec;
intgd2 = Dealias2GLL2*massd2;

%% Spectral Nx to nodal Nxd2
pht = zeros(Nxd2+1,Nx+1);
nx = Nx+1;
lx = Nxd2+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd2(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
Nxspec_to_Nxd2nodal = pht;              % Nx spectral to Nxd2 nodal

GLL2Dealias2 = Nxspec_to_Nxd2nodal*LegendreTransform;

gradm1d2 = GLL2Dealias2*gradm1;

%-------------------------------------------------- 
% Nxd3
pht = zeros(Nx+1,Nxd3+1);
nx = Nxd3+1;
lx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=x(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
dealiased_intpx3 = pht;             % Nxd3 spectral to Nx nodal
%% Build spectral transform matrix on Nxd3 
pht = zeros(Nxd3+1,Nxd3+1);
nx = Nxd3+1;
lx = Nxd3+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd3(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
Nxd3_spec2nodal = pht;              % Nxd spectral to Nxd nodal
Nxd3_nodal2spec = inv(pht);         % Nxd nodal to Nxd spectral

Dealias2GLL3 = dealiased_intpx3*Nxd3_nodal2spec;

intgd3 = Dealias2GLL3*massd3;

%% Spectral Nx to nodal Nxd3
pht = zeros(Nxd3+1,Nx+1);
nx = Nx+1;
lx = Nxd3+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd3(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end
Nxspec_to_Nxd3nodal = pht;         % Nx spectral to Nxd3 nodal

GLL2Dealias3 = Nxspec_to_Nxd3nodal*LegendreTransform;

gradm1d3 = GLL2Dealias3*gradm1;

toc

return



