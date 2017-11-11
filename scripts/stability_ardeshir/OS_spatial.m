function [A,B] = OS_spat(U, D, omega, beta, Re)
%
% [A,B]=OS_temp(U,D,ALPHA,BETA,RE)
% OS_temp sets up the operator corrsponding to local temporal 1D stability
% (Orr-Sommerfeldt) equation.
%
% INPUT
%   U: structure array containing meanflow quantities U.u, U.uy, U.uyy 
%   D: spectral ifferential operator
%   alpha: stremwise wavenumber
%   beta: spanwise wavenumber
%   Re: Reynolds number
%
% OUTPUT
%   A and B: matrices corresponding to the eigenvalue problem A.q = omega.B.q
%
% (c) Ardeshir Hanifi 2014
%


N=length(U.u);
D1=D(:,:,1);
D2=D(:,:,2);
D3=D(:,:,3);
D4=D(:,:,4);
%------------------------------------------------------------------
% Set up operator matrices

zi = sqrt(-1);
k2 = alpha*alpha + beta*beta;
I = eye(N);

%A=(D4-2*k2*D2)*zi/Re+diag(alpha*U.u)*(D2-k2*I)-diag(alpha*U.uyy);
%B=D2-k2*I;

% Zeroth order terms (OS)
A0_V = 1/Re*D4 + (beta^4)/Re - 2*(beta^2)/Re*D2 - zi*omega*(beta^2) - zi*omega*D2;
% First order terms (OS)
A1_V = -4/Re*D3 + 4*(beta^2)/Re - 2*zi*omega*D1 - zi*diag(U.u)*D2 + zi*(beta^2)*diag(U.u) + zi*diag(U.uyy);
% Second order terms (OS)
A2_V = 4/Re*D2 + 2*zi*diag(U.u)*D2;

A0_ETA = zeros(size(A0_V));
A1_ETA = zeros(size(A0_V));
A2_ETA = zeros(size(A0_V));


%% Squire's Equations (Not Used)
% Zeroth order terms (SQ)
B0_V = -zi*beta*diag(U.uy);
B1_V = zeros(size(B0_V));
B2_V = zeros(size(B0_V));

% Zeroth order terms (SQ)
B0_ETA = zi*omega + D2/Re - (beta^2)/Re;
% First order terms (SQ)
B1_ETA = -2/Re*D - zi;
% Second order terms (SQ)
B2_ETA = zeros(size(B0_ETA));


%------------------------------------------------------------------
% Boundary conditions
% Replace the empty rows of B with a multiple of same rows in A. This
% removes the singularity of the system and introduce eigenvalues which are
% located far from the physical ones by approperiate choice of multiple,
eps = 1e-4*zi;
Ov = zeros(1,N);

% v(ymax)=0
A0_V(1,:) = Ov;
A0_V(1,1) = 0;

A1_V(1,:) = Ov;
A1_V(1,1) = 1;

A2_V(1,:) = Ov;
A2_V(1,1) = 1;

% v(0)=0
A0_V(N,:) = Ov;
A0_V(N,N) = 0;

A1_V(N,:) = Ov;
A1_V(N,N) = 1;

A2_V(N,:) = Ov;
A2_V(N,N) = 1;


% dv/dy(ymax)=0
A0_V(2,:) = 1/3*D1(1,:);
A0_V(2,1) = 0;

A1_V(2,:) = 1/3*D1(1,:);

A2_V(2,:) = 1/3*D1(1,:);
A2_V(2,1) = 1;


%---------------------------------------------------------------------- 
% % v(ymax)=0
% A(1,:) = Ov;
% A(1,1) = 1;
% B(1,:) = A(1,:)*eps; 
% 
% % dv/dy(ymax)=0
% A(2,:) = D1(1,:);
% B(2,:) = A(2,:)*eps;
% 
% % dv/dy(0)=0
% A(N-1,:) = D1(N,:);
% B(N-1,:) = A(N-1,:)*eps;
% 
% % v(0) = 0
% A(N,:) = Ov;
% A(N,N) = 1;
% B(N,:) = A(N,:)*eps;

end
