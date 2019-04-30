function [y1, vec, A, B] = PitchingResolvent(uprof,yprof,x0,y0,snx,sny,N,Re,alpha,beta,omega,eta0)

%ifmap=0;
%delta_alpha = 0.2*alpha;
%nalpha=10;
%alpha_range = linspace(alpha-delta_alpha,alpha+delta_alpha,nalpha);

yprof = yprof - min(yprof);
ymax = max(yprof);

% Check stability properties of the profile
[D0,D1,D2,D3,D4] = ChebMat2(N);
Jac = ymax/2;

D1 = (1/Jac)*D1;
D2 = ((1/Jac)^2)*D2;
D3 = ((1/Jac)^3)*D3;
D4 = ((1/Jac)^4)*D4;
phy2cheb = inv(D0);

%% Select the collocation nodes and set up derivative operator
% Highest order of derivative required
M=2;
[y, D] = chebdif(N+1, M);

% y goes from [1, -1]
% y(1) is far-field
% y(N+1) is the wall

y1 = (y+1)*Jac;
u_interp = interp1(yprof,uprof,y1,'pchip');

ucheb = phy2cheb*u_interp;
du = D1*ucheb;
ddu = D2*ucheb;

U.u=u_interp;
U.uy=du;
U.uyy=ddu;

U.w=zeros(size(y)); % This field is required as the code is written for a 3D flow
U.wy=zeros(size(y)); % This field is required as the code is written for a 3D flow

D(:,:,1) = D(:,:,1)*(1/Jac);
D(:,:,2) = D(:,:,2)*((1/Jac)^2);

%% Set up the eigenvalue problem
%dbstop in StabilityProperties at 41
[A,B] = LNS_resolvent(U, D, alpha, beta, Re);

% Generate RHS forcing
[r c] = size(A);
N1=N+1;
% first n rows are the divergence equations.
% First point is the wall point.
%eta0 = 1.;
eta = eta0; %exp(iwt)
zi = sqrt(-1.);
etav = zi*omega*eta0;

Rot = [0 -eta;
       eta 0];

dr = Rot*[x0; y0];
dx = dr(1);
dy = dr(2);
dvx = -etav*y0;
dvy = etav*x0;

% Wall parallel direction
snxp = sny;
snyp = -snx;
if (snxp<0)
  snxp=abs(snxp); % we want it along positive x
  snyp=-snyp;
end  

% Project parallel and normal to the wall
dypp = dx*snx + dy*sny;
dxpp = dx*snxp + dy*snyp;

dvypp = dvx*snx + dvy*sny;
dvxpp = dvx*snxp + dvy*snyp;

dUdy = du(end); 
ux_w = dvxpp - dypp*dUdy; 
uy_w = dvypp;

% RHS forcing
FF = zeros(r,1);
FF(2*N1)=ux_w;
FF(3*N1)=uy_w;

Resolvent = (omega*B - A);
vec = Resolvent\FF;


%[eigvec, eigval] = eig(A,B);
%dbstop in StabilityProperties at 46
%eigval=diag(eigval);


