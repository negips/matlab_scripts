function [eigvec eigval] = StabilityProperties(uprof,yprof,N,Re,alpha,beta,ifOS)

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
if ifOS
  M=4;
else
  M=2;
end
[y, D] = chebdif(N+1, M);

y1 = (y+1)*ymax/2;
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
if (ifOS)
  D(:,:,3) = D(:,:,3)*((1/Jac)^3);
  D(:,:,4) = D(:,:,4)*((1/Jac)^4);
end

%% Set up the eigenvalue problem
%dbstop in StabilityProperties at 41
if (ifOS)
  [A,B] = OS_temp(U, D, alpha, beta, Re);
else
  [A,B] = LNS_temp(U, D, alpha, beta, Re);
end
%% Solve the eigenvalue problem
[eigvec, eigval] = eig(A,B);
%dbstop in StabilityProperties at 46
eigval=diag(eigval);


