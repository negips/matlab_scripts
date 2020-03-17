% This script runs the temporal 1D stability problem
%
% (c) Ardeshir Hanifi 2014
%
clear all
clc
close all
%% Numerical parameters
N = 128; % number of Cehebychev modes

%% Flow parameters
%Re=54385;
%alpha=0.1555;
%beta=0.0;

Re=1000;
alpha=1.0;
beta=1.0;

%% Select the collocation nodes and set up derivative operator
M=2; % Highest order of derivative required
[y, D] = chebdif(N, 2);

%% Velocity profile corresponding to Plane Poiseuille Flow
U.u=1-y.^2;
U.uy=-2*y;
U.w=zeros(size(y)); % This field is required as the code is written for a 3D flow
U.wy=zeros(size(y)); % This field is required as the code is written for a 3D flow

ymax = 2;
%y1 = (y+1)*ymax/2;
%U.u=1-exp(V0*Re*y1);
%U.uy=-V0*Re*exp(V0*Re*y1);

D(:,:,1) = D(:,:,1)/(ymax)*2;
D(:,:,2) = D(:,:,2)/(ymax^2)*(2^2);

%% Set up the eigenvalue problem
[A,B] = LNS_temp(U, D, alpha, beta, Re);

%% Solve the eigenvalue problem
[eigvec, eigval] = eig(A,B);
eigval=diag(eigval);

%% Plot spectrum and eigenvectors
header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
              {', beta= '},{num2str(beta)});

figure(1)
plot(real(eigval),imag(eigval),'o');
title(header);
ylabel('$\omega_{i}$')
xlabel('$\omega_{r}$')

%plot_LNS(eigval,eigvec,y,header)
%xlim([0 0.2]);
%ylim([-0.2 0])
