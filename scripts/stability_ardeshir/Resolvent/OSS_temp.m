function [A,B] = OSS_temp(U, D, alpha, beta, Re)
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
D4=D(:,:,4);
%------------------------------------------------------------------
% Set up operator matrices

zi = sqrt(-1);
k2 = alpha*alpha + beta*beta;
I = eye(N);
z0 =zeros(N);

A11=(D4-2*k2*D2+k2*k2*I)*zi/Re+diag(alpha*U.u)*(D2-k2*I)-diag(alpha*U.uyy);
A21=beta*diag(U.uy);
A22=(D2-k2*I)*zi/Re+diag(alpha*U.u);

B11=D2-k2*I;

A=[A11 z0; A21 A22];
B=[B11 z0;z0 I];
%------------------------------------------------------------------
% Boundary conditions
% Replace the empty rows of B with a multiple of same rows in A. This
% rmoves the singularity of the system and introduce eigenvalues which are
% located far from the physical ones by approperiate choice of multiple,
eps = 1e-4*zi;
Ov = zeros(1,N*2);

% v(ymax)=0
A(1,:) = Ov;
A(1,1) = 1;
B(1,:) = A(1,:)*eps; 

% dv/dy(ymax)=0
A(2,1:N) = D1(1,:);
B(2,:) = A(2,:)*eps;

% dv/dy(0)=0
A(N-1,1:N) = D1(N,:);
B(N-1,:) = A(N-1,:)*eps;

% v(0) = 0
A(N,:) = Ov;
A(N,N) = 1;
B(N,:) = A(N,:)*eps;

% eta(ymax)=0
A(N+1,:)   = Ov;
A(N+1,N+1) = 1;
B(N+1,:)   = A(N+1,:)*eps; 

% eta(0) = 0
A(2*N,:)   = Ov;
A(2*N,2*N) = 1;
B(2*N,:)   = A(2*N,:)*eps;

end