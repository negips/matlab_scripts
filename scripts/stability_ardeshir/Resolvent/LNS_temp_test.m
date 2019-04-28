function [A,B] = LNS_temp(U, D, alpha, beta, Re)
%
% LNS_temp sets up the operator of general eigenvalue problem A.q = omega.B.q
% corrsponding to local temporal 1D stability equations.
%
% [A,B]=LNS_tmp(U,D,ALPHA,BETA,RE).
% 
% Input:
%   U:  structure containig streamwise and spanwise velocities and their 
%       first normal derivatives U.u, U.uy, U.w, U,wy
%   D:  sepctral differential operator
%   ALPHA: streamwise wavenumber
%   BETA: spanwise wavenumber
%   RE: Reynolds number
%
% OUTPUT
%   A and B: matrices corresponding to the eigenvalue problem A.q = omega.B.q
%
% (c) Ardeshir Hanifi, 2014

N=length(U.u);
D1=D(:,:,1);
D2=D(:,:,2);
%------------------------------------------------------------------
% Set up operator matrices

zi = sqrt(-1);
I = eye(N);
O = zeros(N);
xi=(D2-(beta^2+alpha^2)*I)/Re - zi*diag(alpha*U.u + beta*U.w);

A11 = zi*alpha*I;
A12 = D1;
A13 = zi*beta*I;
A21 = xi;
A22 = -diag(U.uy);
A24 = -zi*alpha*I;
A32 = xi;
A34 = -D1;
A42 = -diag(U.wy);
A43 = xi;
A44 = -zi*beta*I;

A = [A11, A12, A13, O  ; ...
     A21, A22, O,   A24; ...
     O,   A32, O,   A34; ...
     O,   A42, A43, A44];

%------------------------------------------------------------------

B = [O,     O,   O,     O; ...
     -zi*I, O,   O,     O; ...
     O,   -zi*I, O,     O; ...
     O,     O,   -zi*I, O];

%------------------------------------------------------------------
% Boundary conditions

Ov = zeros(1,4*N);
eps=1i*1e-8*0;
% u
A(N+1,:) = Ov;
A(N+1,1) = 1;
B(N+1,:) = A(N+1,:)*eps; 

A(2*N,:) = Ov;
A(2*N, N) = 1;
B(2*N,:) = A(2*N,:)*eps;

%v
A(2*N+1,:) = Ov;
A(2*N+1, N+1) = 1;
B(2*N+1,:) = A(2*N+1,:)*eps;

A(3*N,:) = Ov;
A(3*N, 2*N) = 1;
B(3*N,:) = A(3*N,:)*eps;

%w
A(3*N+1,:) = Ov;
A(3*N+1, 2*N+1) = 1;
B(3*N+1,:) = A(3*N+1,:)*eps;

A(4*N,:) = Ov;
A(4*N, 3*N) = 1;
B(4*N,:) = A(4*N,:)*eps;

% %p
% A(1,:) = Ov;
% A(1,N+1:2*N) = D2(1,:)/Re;
% A(1,3*N+1:4*N) = -D1(1,:);
% B(1,:) = A(1,:)*eps;
% 
% A(N,:) = Ov;
% A(N,N+1:2*N) = D2(N,:)/Re;
% A(N,3*N+1:4*N) = -D1(N,:);
% B(N,:) = A(N,:)*eps;
end