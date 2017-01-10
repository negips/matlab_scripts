% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% dynamic fault with slip weakening,
% paraxial absorbing boundary conditions
% and zero initial conditions
% in a structured undeformed grid.
%
% Version 2: domain = rectangular
%            medium = general (heterogeneous)
%            boundaries = 1 fault + 3 paraxial
%	     time scheme = Newmark as in SPECFEM3D: alpha=1,beta=0,gamma=1/2
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here the parameters of the square box domain and mesh : ****
LX=50e3;
LY=50e3/3;
%NELX = 75; NELY = 25; P = 8; % polynomial degree
NELX = 150; NELY = 50; P = 4; % polynomial degree
%********

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
x = x-LX/2;
nglob = length(x);

RHO = 2670.;
VS  = 3464.;

%------------------------------------------
% STEP 2: INITIALIZATION

[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

W     = zeros(NGLL,NGLL,NEL);	% for internal forces
M     = zeros(nglob,1);		% global mass matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored

%**** Set here the parameters of the time solver : ****
NT = 2200; % number of timesteps
CFL   = 0.6; 			% stability number = CFL_1D / sqrt(2)
%********

dt    = Inf;  			% timestep (set later)

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi^2 ;
coefint2 = jac/dy_deta^2 ;

% FOR EACH ELEMENT ...
for ey=1:NELY, 
for ex=1:NELX, 

  e = (ey-1)*NELX+ex;
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
  rho(:,:) = RHO;
  mu(:,:)  = RHO* VS^2;
%********

 % Diagonal mass matrix
  M(ig) = M(ig) + wgll2 .*rho *jac;

 % Local contributions to the stiffness matrix K
 %  WX(:,:,e) = wgll2 .* mu *jac/dx_dxi^2;
 %  WY(:,:,e) = wgll2 .* mu *jac/dy_deta^2;
  W(:,:,e) = wgll2 .* mu;

 % The timestep dt is set by the stability condition
 %   dt = CFL*min(dx/vs)
  vs = sqrt(mu./rho); 
  if dxe<dye
    vs = max( vs(1:NGLL-1,:), vs(2:NGLL,:) );
    dx = repmat( diff(xgll)*0.5*dxe ,1,NGLL); 
  else
    vs = max( vs(:,1:NGLL-1), vs(:,2:NGLL) );
    dx = repmat( diff(xgll)'*0.5*dye ,NGLL,1); 
  end
  dtloc = dx./vs;
  dt = min( [dt dtloc(1:end)] );

end
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;


%-- Absorbing boundaries (first order): 
impedance = RHO*VS;
% Left
ng = NELY*(NGLL-1)+1;
BcLeftIglob = zeros(ng,1);
BcLeftC = zeros(ng,1);
for ey=1:NELY,
  ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
  e=(ey-1)*NELX+1;
  BcLeftIglob(ip) = iglob(1,1:NGLL,e);
  BcLeftC(ip) = BcLeftC(ip) + dy_deta*wgll*impedance ;
end
% Right
ng = NELY*(NGLL-1)+1;
BcRightIglob = zeros(ng,1);
BcRightC = zeros(ng,1);
for ey=1:NELY,
  ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
  e=(ey-1)*NELX+NELX;
  BcRightIglob(ip) = iglob(NGLL,1:NGLL,e);
  BcRightC(ip) = BcRightC(ip) + dy_deta*wgll*impedance ;
end
% Top
ng = NELX*(NGLL-1)+1;
BcTopIglob = zeros(ng,1);
BcTopC = zeros(ng,1);
for ex=1:NELX,
  ip = (NGLL-1)*(ex-1)+[1:NGLL] ;
  e=(NELY-1)*NELX+ex;
  BcTopIglob(ip) = iglob(1:NGLL,NGLL,e);
  BcTopC(ip) = BcTopC(ip) + dx_dxi*wgll*impedance ;
end

% The mass matrix needs to be modified at the boundary
% for the IMPLICIT treatment of the term C*v.
% Fortunately C is diagonal.
M(BcLeftIglob)  = M(BcLeftIglob)  +half_dt*BcLeftC;
M(BcRightIglob) = M(BcRightIglob) +half_dt*BcRightC;
M(BcTopIglob)   = M(BcTopIglob)   +half_dt*BcTopC;


%-- DYNAMIC FAULT at bottom boundary
FaultNglob = NELX*(NGLL-1)+1;
FaultIglob = zeros(FaultNglob, 1); 
FaultB = zeros(FaultNglob, 1); 
for ex=1:NELX,
  ip = (NGLL-1)*(ex-1)+[1:NGLL];
  e = ex;
  FaultIglob(ip) = iglob(1:NGLL,1,e);
  FaultB(ip) = FaultB(ip) + dx_dxi*wgll;
end
FaultZ = M(FaultIglob)./FaultB /half_dt;
FaultX = x(FaultIglob);
FaultV = zeros(FaultNglob,NT);
FaultD = zeros(FaultNglob,NT);
% background
FaultNormalStress = 120e6;
FaultInitStress = repmat(70e6,FaultNglob,1);
FaultState = zeros(FaultNglob,1);
FaultFriction.MUs = repmat(0.677,FaultNglob,1);
FaultFriction.MUd = repmat(0.525,FaultNglob,1);
FaultFriction.Dc  = 0.4;
% barrier
isel = find(abs(FaultX)>15e3);
FaultFriction.MUs(isel) = 1e4; % barrier
FaultFriction.MUd(isel) = 1e4; % barrier
% nucleation
isel = find(abs(FaultX)<=1.5e3);
FaultInitStress(isel) = 81.6e6;
FaultFriction.W = (FaultFriction.MUs-FaultFriction.MUd)./FaultFriction.Dc;
FaultStrength = friction(FaultState,FaultFriction)*FaultNormalStress ...
                - FaultInitStress; % strength excess


%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTxseis = [-16e3:600:16e3]';		% x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
OUTyseis = repmat(7.5e3,OUTnseis,1);	% y coord of receivers
%********
% receivers are relocated to the nearest node
% OUTdseis = distance between requested and relocated receivers
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 50;
OUTit = 0;
OUTindx = Init2dSnapshot(iglob);



%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark scheme with
% alpha=1, beta=0, gamma=1/2
%
for it=1:NT,

 % update
  d = d + dt*v + half_dt_sq*a; 
 % prediction 
  v = v + half_dt*a;
  a(:) = 0; 

 % internal forces -K*d(t+1) 
 % stored in global array 'a'
  for e=1:NEL,
   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d(ig);
   %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;	
    d_eta = local*H;
   %element contribution to internal forces
   %local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
    wloc = W(:,:,e);
    d_xi = wloc.*d_xi;
    d_xi = H * d_xi;
    d_eta = wloc.*d_eta;
    d_eta = d_eta *Ht;
    local = coefint1* d_xi  + coefint2* d_eta ;
   %assemble into global vector
    a(ig) = a(ig) -local;
  end 

 % absorbing boundaries:
  a(BcLeftIglob)  = a(BcLeftIglob)  - BcLeftC  .* v(BcLeftIglob);
  a(BcRightIglob) = a(BcRightIglob) - BcRightC .* v(BcRightIglob);
  a(BcTopIglob)   = a(BcTopIglob)   - BcTopC   .* v(BcTopIglob) ;

 % fault boundary condition: slip weakening
  FaultVFree = v(FaultIglob) + half_dt*a(FaultIglob)./M(FaultIglob);
  TauStick = FaultZ .*FaultVFree;
 % TauStick = a(FaultIglob)./FaultB;
  Tau = min(TauStick,FaultStrength); 
  a(FaultIglob) = a(FaultIglob) - FaultB .*Tau;

 % second pass in EXPLICIT implementation of absorbing boundary conditions
 % (does not require modification of M)
 % BUT DOES NOT WORK !! (nearly as unstable as explicit with single pass)
%  a(BcLeftIglob) = a(BcLeftIglob) - half_dt*BcLeftC.*a(BcLeftIglob)./M(BcLeftIglob);
%  a(BcRightIglob) = a(BcRightIglob) - half_dt*BcRightC.*a(BcRightIglob)./M(BcRightIglob);
%  a(BcTopIglob) = a(BcTopIglob) - half_dt*BcTopC.*a(BcTopIglob)./M(BcTopIglob);

 % solve for a_new:
  a = a ./M ;

 % correction
  v = v + half_dt*a;

  FaultState = max(2*d(FaultIglob),FaultState);
  FaultStrength = friction(FaultState,FaultFriction)*FaultNormalStress ...
                  -FaultInitStress;
%% rupture nucleation through time weakening, like Andrews 76
%  VNUC = 0.4;
%%  if time(it)<10/VNUC
%    %FaultStrength = min( FaultStrength,...
%    ix = find(abs(FaultX)<=10);
%    FaultStrength(ix) = min(FaultStrength(ix), ...
%              max(0.5,0.55+(abs(FaultX(ix))-VNUC*time(it))*0.05/(1.0*dxe))...
%                     - FaultInitStress(ix) ) ;
%%  end
  FaultV(:,it) = 2*v(FaultIglob);
  FaultD(:,it) = 2*d(FaultIglob);

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[0 2]);
    hold on
    plot(OUTxseis,OUTyseis,'^')
    hold off

    drawnow

  end

end % ... of time loop

figure(3)
sdt=10;
sdx=4;
surf(time(1:sdt:end),FaultX(1:sdx:end),FaultV(1:sdx:end,1:sdt:end))
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')
