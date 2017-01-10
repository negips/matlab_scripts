% SEM1D	applies the Spectral Element Method
% to solve the 1D SH wave equation, 
% with stress free boundary conditions,
% zero initial conditions
% and a time dependent force source.
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%
% THIS VERSION: implicit generalized alpha scheme
%		(J. Chung and G.M. Hulbert, JAM 1993)
%		which allows for a user-controlled dissipation at high-frequency 
%		while mantaining low dissipation at low-frequency.
%		As it is an implicit scheme, large timesteps are allowed.
%		However, it introduces dispersion. 
%		In this version the inversion of the dynamic matrix 
%		is done by brute force.

%------------------------------------------
% STEP 1: MESH GENERATION
% The interval [0,L] is divided into NEL non overlapping elements
% The elements are defined by their "control nodes" X
L=20;
NEL = 20;
X = [0:NEL]'*L/NEL;

%------------------------------------------
% STEP 2: INITIALIZATION

P = 6; % polynomial degree
NGLL = P+1; % number of GLL nodes per element
NT = 125; % number of timesteps

ABSO_TOP = 0;
ABSO_BOTTOM =0;

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL);

iglob = zeros(NGLL,NEL);	% local to global index mapping
rho = zeros(NGLL,NEL);		% density
mu = zeros(NGLL,NEL);		% shear modulus
nglob = NEL*(NGLL-1) + 1;	% number of global nodes
coor = zeros(nglob,1);		% coordinates of GLL nodes
M = zeros(nglob,1);		% global mass matrix, diagonal
K = zeros(nglob,nglob);		% global stiffness matrix
Kloc = zeros(NGLL,NGLL,NEL);	% local stiffness matrix
CFL = 4; 			% stability number, can be high (but introduces dispersion)
dt = Inf;  			% timestep (set later)

r = 1.; 	% dissipation parameter of the generalized-alpha scheme
		% 1 = no dissipation
		% 0 = maximal high-frequency dissipation
alpha_m = (2*r-1)/(r+1);
alpha_f = r/(1+r);
beta = 1/4*(1-alpha_m+alpha_f)^2;
gamma = 0.5-alpha_m+alpha_f;

for e=1:NEL, % FOR EACH ELEMENT ...

 % The table I = iglob(i,e) maps the local numbering of the 
 % computational nodes (i,e) to their global numbering I.
 % 'iglob' is used to build global data 
 % from local data (contributions from each element)
  iglob(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';

 % Coordinates of the computational (GLL) nodes
  dxe = X(e+1)-X(e);
  coor(iglob(:,e)) = 0.5*(X(e)+X(e+1)) + 0.5*dxe*xgll ;

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
  rho(:,e) = 1;
  mu(:,e) = 1;

 % For this simple mesh the jacobian of the global-local
 % coordinate map is a constant
  dx_dxi = 0.5*dxe;

 % The diagonal mass matrix is stored in a global array. 
 % It is assembled here, from its local contributions
 % Nodes at the boundary between two elements get 
 % contributions from both.
  ig = iglob(:,e);
  M(ig) = M(ig) + wgll .*rho(:,e) *dx_dxi;

% The stiffness matrix K is not assembled at this point
% (it is sparse, block-diagonal)
% We only store its local contributions
  W = mu(:,e).*wgll/dx_dxi;
%  K(:,:,e) = H * diag(W)* H';
  Kloc(:,:,e) = H * ( repmat(W,1,NGLL).* H');
  K(ig,ig) = K(ig,ig) + Kloc(:,:,e);

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
  vs = sqrt(mu(:,e)./rho(:,e)); 
  vs = max( vs(1:NGLL-1), vs(2:NGLL) );
  dx = abs(diff( coor(iglob(:,e)) ));
  dt = min(dt, min(dx./vs));
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;

D = (1-alpha_m)*diag(M) + ((1-alpha_f)*beta*dt^2)*K ;
%invD = inv( (1-alpha_m)*diag(M) + ((1-alpha_f)*beta*dt^2)*K );

% absorbing boundaries
% The mass matrix needs to be modified at the boundary
% for the implicit treatment of the term C*v.
if ABSO_TOP
  BcTopC = sqrt(rho(NGLL,NEL)*mu(NGLL,NEL));
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
  M(1) = M(1)+half_dt*BcBottomC;
end

% Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);
lhs = zeros(nglob,1);
dlhs = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
% or velocity amplitude of an incoming wave
F_IS_WAVE = 0;
if ~F_IS_WAVE,
  Fx = L/2;
  [Fdist,Fix] = min( abs(coor-Fx) );
end
Ff0 = 0.25; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = ricker( (1:NT)'*dt-0.5*dt, Ff0,Ft0);

% output arrays
OUTdt = 1;	% output every OUTdt timesteps
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output receivers at these locations
OUTx  = [0:dxe:L]';
OUTnx = length(OUTx);
OUTix = zeros(OUTnx,1);
% relocate to nearest GLL node
OUTdist = zeros(OUTnx,1);
for i=1:OUTnx,
  [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) );
end
OUTx = coor(OUTix);
OUTd = zeros(OUTnx,OUTnt);
OUTv = zeros(OUTnx,OUTnt);
OUTa = zeros(OUTnx,OUTnt);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
%

for it=1:NT,

  dlhs = d;
  d = d +dt*v + (dt^2*(0.5-beta))*a;
  dlhs = (1-alpha_f)*d + alpha_f*dlhs;
  v = v +(dt*(1-gamma))*a;

 % internal forces at mid-step -K*d(t+1/2) 
  lhs(:) = 0; % store -K*d in a global array
  for e=1:NEL,
    ix = iglob(:,e);
    lhs(ix) = lhs(ix) - Kloc(:,:,e)*dlhs(ix) ;
  end 
  lhs = lhs - alpha_m*M.*a;

 % add external forces
  if ~F_IS_WAVE, lhs(Fix) = lhs(Fix) + Ft(it); end

% % absorbing boundary
%  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
%  if ABSO_BOTTOM, 
%    if F_IS_WAVE % incident wave, from bottom
%     %boundary term = -impedance*v_outgoing + impedance*v_incoming
%      a(1) = a(1) -BcBottomC*( v(1) -Ft(it) ) + BcBottomC*Ft(it);
%    else
%      a(1) = a(1) -BcBottomC*v(1);
%    end
%  end

 % acceleration: a = (-K*d +F)/M
  %a = invD*lhs ;
  a = D\lhs ;

 % update
  v = v + (gamma*dt)*a;
  d = d + (beta*dt^2)*a;


%------------------------------------------
% STEP 4: OUTPUT
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;
    OUTd(:,OUTit) = d(OUTix);
    OUTv(:,OUTit) = v(OUTix);
    OUTa(:,OUTit) = a(OUTix);
  end

end % of time loop

t = (1:OUTit) *OUTdt*dt;
PlotSeisTrace(OUTx,t,OUTv);
