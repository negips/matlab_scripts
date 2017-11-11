function [A,B] = OS_temp(U, D, alpha, beta, Re)
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

A=(D4-2*k2*D2)*zi/Re+diag(alpha*U.u)*(D2-k2*I)-diag(alpha*U.uyy);
B=D2-k2*I;
%------------------------------------------------------------------
% Boundary conditions
% Replace the empty rows of B with a multiple of same rows in A. This
% rmoves the singularity of the system and introduce eigenvalues which are
% located far from the physical ones by approperiate choice of multiple,
eps = 1e-4*zi;
Ov = zeros(1,N);

% v(ymax)=0
A(1,:) = Ov;
A(1,1) = 1;
B(1,:) = A(1,:)*eps; 

% dv/dy(ymax)=0
A(2,:) = D1(1,:);
B(2,:) = A(2,:)*eps;

% dv/dy(0)=0
A(N-1,:) = D1(N,:);
B(N-1,:) = A(N-1,:)*eps;

% v(0) = 0
A(N,:) = Ov;
A(N,N) = 1;
B(N,:) = A(N,:)*eps;

end
