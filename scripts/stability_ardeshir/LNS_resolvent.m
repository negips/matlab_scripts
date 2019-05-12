function [A,B] = LNS_resolvent(U, D, alpha, beta, Re)
%
% LNS_resolvent sets up the operator of general eigenvalue problem A.q = omega.B.q
% corrsponding to local temporal 1D stability equations, with Identity at the boundary conditions so that the forcing dictates this boundary value
%
% [A,B]=LNS_resolvent(U,D,ALPHA,BETA,RE).
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
%Iw = I*omega;

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

% u
% far-field Dirichlet
%A(N+1,:) = Ov;
%A(N+1,1) = 1;
%B(N+1,:) = Ov;

% far-field Neumann
A(N+1,:) = Ov;
A(N+1,1:N) = D1(1,:);
%A(N+1,1) = 1;
B(N+1,:) = Ov;

% wall Dirichlet
A(2*N,:) = Ov;
A(2*N, N) = 1;
B(2*N,:) = Ov;

%v
% far-field Dirichlet
%A(2*N+1,:) = Ov;
%A(2*N+1, N+1) = 1;
%B(2*N+1,:) = Ov;

% far-field Neumann
A(2*N+1,:) = Ov;
A(2*N+1,N+1:2*N) = D1(1,:);
%A(2*N+1, N+1) = 1;
B(2*N+1,:) = Ov;


% wall Dirichlet
A(3*N,:) = Ov;
A(3*N, 2*N) = 1;
B(3*N,:) = Ov;

%w
% far-field Dirichlet
A(3*N+1,:) = Ov;
A(3*N+1, 2*N+1) = 1;
B(3*N+1,:) = Ov;

% wall Dirichlet
A(4*N,:) = Ov;
A(4*N, 3*N) = 1;
B(4*N,:) = Ov;

% Remove rows,columns
%dels = [2*N 3*N 4*N];
dels = [];

A(dels,:)=[];
A(:,dels)=[]; 

B(dels,:)=[];
B(:,dels)=[]; 

end





