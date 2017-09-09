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

% THIS VERSION: applied force on the bottom boundary
%		absorbing condition on the top boundary
% => analyze numerical impedance

%------------------------------------------
% STEP 1: MESH GENERATION
% The interval [0,L] is divided into NEL non overlapping elements
% The elements are defined by their "control nodes" X
L=10;
NEL = 10;
X = [0:NEL]'*L/NEL;

% add a viscous layer (one element) next to the boundary :
NEL_VISC = 0;
tvisc = 0.02; % viscous time /dt

%------------------------------------------
% STEP 2: INITIALIZATION

P = 8; % polynomial degree
NGLL = P+1; % number of GLL nodes per element
NT = 4096; % number of timesteps

ABSO_TOP = 1;
ABSO_BOTTOM =0;

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL,'gll');

iglob = zeros(NGLL,NEL);	% local to global index mapping
rho = zeros(NGLL,NEL);		% density
mu = zeros(NGLL,NEL);		% shear modulus
nglob = NEL*(NGLL-1) + 1;	% number of global nodes
coor = zeros(nglob,1);		% coordinates of GLL nodes
M = zeros(nglob,1);		% global mass matrix, diagonal
K = zeros(NGLL,NGLL,NEL);	% local stiffness matrix
CFL = 0.5; 			% stability number 0.85
dt = Inf;  			% timestep (set later)

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
  M(iglob(:,e)) = M(iglob(:,e)) + wgll .*rho(:,e) *dx_dxi;

% The stiffness matrix K is not assembled at this point
% (it is sparse, block-diagonal)
% We only store its local contributions
  W = mu(:,e).*wgll/dx_dxi;
%  K(:,:,e) = H * diag(W)* H';
  K(:,:,e) = H * ( repmat(W,1,NGLL).* H');

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
  vs = sqrt(mu(:,e)./rho(:,e)); 
  vs = max( vs(1:NGLL-1), vs(2:NGLL) );
  dx = abs(diff( coor(iglob(:,e)) ));
  dt = min(dt, min(dx./vs));
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;

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

% treat implicitly the diagonal part of the damping matrix
%if NEL_VISC
%  Kdiag = zeros(NGLL,NEL_VISC);
%  for e=1:NEL_VISC,
%    ix = iglob(:,e);
%    Kdiag(:,e) = diag(K(:,:,e));
%    M(ix) = M(ix) + 0.1* 0.5*tvisc*Kdiag(:,e);
%  end
%end

% Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);
f = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
% or velocity amplitude of an incoming wave
F_IS_WAVE = 0;
if ~F_IS_WAVE,
  Fx = 0;
  [Fdist,Fix] = min( abs(coor-Fx) );
end
Ff0 = 3; % fundamental frequency
Ft0 = 2/Ff0; % delay
% source time function (at mid-steps)
%Ft = ricker( (1:NT)'*dt-0.5*dt, Ff0,Ft0);
%Ft = gabor( (1:NT)'*dt-0.5*dt, Ff0,Ft0, 80 );
Ft = gaussian( (1:NT)'*dt-0.5*dt, Ff0,Ft0 );

% output arrays
OUTdt = 1;	% output every OUTdt timesteps
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output receivers at these locations
OUTx  = [0:dxe/2:L]';
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
% Explicit Newmark-alpha scheme with
% alpha=1/2, beta=1/2, gamma=1
%

for it=1:NT,

 % prediction of mid-step displacement:
 % d_mid = d_old + 0.5*dt*v_old
  d = d + half_dt*v; 

  f(:) = 0;

 % internal forces at mid-step -K*d(t+1/2) 
  for e=1:NEL,
    ix = iglob(:,e);
    f(ix) = f(ix) - K(:,:,e)*d(ix) ;
  end 

 % viscous layer:
  for e=1:NEL_VISC
    ix = iglob(:,e);
    f(ix) = f(ix) - tvisc* K(:,:,e)*v(ix); 
  end

 % add external forces
  if ~F_IS_WAVE, f(Fix) = f(Fix) + Ft(it); end

 % absorbing boundary
  if ABSO_TOP, f(nglob) = f(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
     %boundary term = -impedance*v_outgoing + impedance*v_incoming
      f(1) = f(1) -BcBottomC*( v(1) -Ft(it) ) + BcBottomC*Ft(it);
    else
      f(1) = f(1) -BcBottomC*v(1);
    end
  end

 % acceleration: a = (-K*d +F)/M
  a = f ./M ;

 % update
 % v_new = v_old + dt*a_new;
 % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
 %       = d_mid + 0.5*dt*v_new
  v = v + dt*a;
  d = d + half_dt*v;


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
figure(1)
PlotSeisTrace(OUTx,t,OUTv);

figure(2)
ista=1;
subplot(311)
[freq,fv]=plot_spec(OUTv(ista,:),dt);
ylabel('velocity bottom boundary')
subplot(312)
[freq,fFt] = plot_spec(Ft,dt);
ylabel('source time function')
hold on; loglog( [freq(1) freq(end)], [1e-3 1e-3]*max(fFt),':');hold off
subplot(313)
loglog(freq,fFt./fv)
ylabel('impedance = f/v')
