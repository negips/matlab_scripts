function [A_os,B_os] = OS_spatial_extended(U, D, omega, beta, Re)
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
I = eye(N);
Z0 = zeros(N);

%A=(D4-2*k2*D2)*zi/Re+diag(alpha*U.u)*(D2-k2*I)-diag(alpha*U.uyy);
%B=D2-k2*I;

%% Expressions from Schmid & Henningson (2000) Transition in Shear Flows. pp 259


% Zeroth order terms (OS)
R0 = zi*omega*D2  - zi*omega*(beta^2)*I + D4/Re - 2*(beta^2)/Re*D2 + (beta^4)/Re;
% First order terms (OS)
R1 = -2*zi*omega*D1 - 4/Re*D3 + 4/Re*(beta^2)*D1 - zi*diag(U.u)*D2 + zi*(beta^2)*diag(U.u) + zi*diag(U.uyy);
% Second order terms (OS)
R2 = 4/Re*D2 + 2*zi*diag(U.u)*D1;

A_os = [-R1  -R0; ...
         I    Z0];

B_os = -[R2  Z0; ...
         Z0  I];


%% Squire's Equations
%% NOT USED/TESTED

S = zi*beta*diag(U.uy);

T0 = -zi*omega*I - D2/Re + (beta^2)/Re;
T1 = 2/Re*D1 - zi*diag(U.u);


A_osq = [-R1  -R0   Z0; ...
          I    Z0   Z0; ...
          Z0  -S   -T0];

B_osq = [R2    Z0   Z0; ...
         Z0    I    Z0; ...
         Z0    Z0   T1];


%------------------------------------------------------------------
% Boundary conditions
% Replace the empty rows of B with a multiple of same rows in A. This
% removes the singularity of the system and introduce eigenvalues which are
% located far from the physical ones by approperiate choice of multiple,
eps = 1e-2;
Ov = zeros(1,N);
Ov2 = zeros(1,2*N);
%Ov3 = zeros(1,3*N);

%---------------------------------------------------------------------- 
% In principle we need 7 boundary conditions.
% Since the equations are 4th order for V and 3rd order for (alpha*V)
% Probably the cause of so many spurious modes

% % V(ymax)=0
A_os(1,:) = Ov2;
A_os(1,N+1) = 1*eps;
B_os(1,:) = A_os(1,:)*eps; 
%B_os(1,:) = Ov2; 

% d(V)/dy(ymax)=0
A_os(2,:) = [Ov D1(1,:)];
B_os(2,:) = A_os(2,:)*eps;
%B_os(2,:) = Ov2;

% V(0) = 0
A_os(N,:) = Ov2;
A_os(N,2*N) = 1*eps;
B_os(N,:) = A_os(N,:)*eps;
%B_os(N,:) = Ov2;

% d(V)/dy(0)=0
A_os(N-1,:) = [Ov D1(N,:)];
B_os(N-1,:) = A_os(N-1,:)*eps;
%B_os(N-1,:) = Ov2;

%% aV
% % aV(ymax)=0
%A_os(3,:) = Ov2;
%A_os(3,1) = 1*eps;
%B_os(3,:) = A_os(1,:)*eps; 
%%B_os(1,:) = Ov2; 
%
%% d(aV)/dy(ymax)=0
%A_os(4,:) = [D1(1,:) Ov];
%B_os(4,:) = A_os(4,:)*eps;
%%B_os(2,:) = Ov2;
%
%% aV(0) = 0
%A_os(N-2,:) = Ov2;
%A_os(N-2,N) = 1*eps;
%B_os(N-2,:) = A_os(N-2,:)*eps;
%%B_os(N,:) = Ov2;
%
%% d(aV)/dy(0)=0
%A_os(N-3,:) = [D1(N,:) Ov];
%B_os(N-3,:) = A_os(N-3,:)*eps;
%B_os(N-3,:) = Ov2;



%% For the bigger Matrix



return
