% Copyright 2016, All Rights Reserved
% Code by Steven L. Brunton
clear all, close all, clc
figpath = './figures/';
addpath('./utils');

% generate Data
sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;
n = 3;
x0=[-8; 8; 27];  % Initial condition

% Integrate
dt = 0.001;
tspan=[dt:dt:200];
N = length(tspan);
%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
%[t,xdat]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);

stackmax = 300;  % the number of shift-stacked rows
lambda = 0;   % threshold for sparse regression (use 0.02 to kill terms)
rmax = 45;       % maximum singular vectors to include

query = importdata('query.out');

nsamples = size(query, 1);
subsamples = 1:25:nsamples;
xdat = query(subsamples,4);
xdat = xdat;
t = query(subsamples,1)*1e-3;
dt = t(2)-t(1);

%% EIGEN-TIME DELAY COORDINATES
clear V, clear dV, clear H
H = zeros(stackmax,size(xdat,1)-stackmax);
for k=1:stackmax
    H(k,:) = xdat(k:end-stackmax-1+k,1);
end
[U,S,V] = svd(H,'econ');
sigs = diag(S);
beta = size(H,1)/size(H,2);
%thresh = optimal_SVHT_coef(beta,0) * median(sigs);
sigall = sum(sigs);
rsum = cumsum(sigs);
%r = find(rsum>0.9*sigall,1);
r = length(sigs);

%r = length(sigs(sigs>thresh));
%r=min(rmax,r);

%figure(10)
%i=5;
%plot3(V(:,i),V(:,i+1),V(:,i+2))

%break

%% COMPUTE DERIVATIVES
% compute derivative using fourth order central difference
% use TVRegDiff if more error 
dV = zeros(length(V)-5,r);
for i=3:length(V)-3
    for k=1:r
        dV(i-2,k) = (1/(12*dt))*(-V(i+2,k)+8*V(i+1,k)-8*V(i-1,k)+V(i-2,k));
    end
end  
% concatenate
x = V(3:end-3,1:r);
dx = dV;
[nx ny] = size(x);
xdat2 = xdat(3:nx+2);

%%  BUILD HAVOK REGRESSION MODEL ON TIME DELAY COORDINATES
% This implementation uses the SINDY code, but least-squares works too
% Build library of nonlinear time series
polyorder = 1;
Theta = poolData(x,r,polyorder,0);

% normalize columns of Theta (required in new time-delay coords)
for k=1:size(Theta,2)
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end 
m = size(Theta,2);
% compute Sparse regression: sequential least squares
% requires different lambda parameters for each column
clear Xi
%for k=1:r-1
%    Xi(:,k) = sparsifyDynamics(Theta,dx(:,k),lambda*k,1);  % lambda = 0 gives better results 
%end

% get projections onto modes
for k=1:r
    [pXi(:,k) resid(:,k)]= CalcProjections(Theta,dx(:,k)); 
end

i=10; figure(1); plot(x(:,i)); figure(2): plot(resid(:,i),'r'); %figure(3); plot(pXi(:,i))

cols = lines(10);

figure(4)
plot(xdat2, 'k'); hold on
for j=1:1
  i=20+j;
  vari = var(x(:,i));
  xmean = mean(x(:,i));
  tol = 0.1*max(x(:,i));    
  ind = find(abs(x(:,i))<tol);
  xdat3=xdat2;    
  xdat3(ind) = nan;
  plot(xdat3, 'o ', 'Color', cols(j,:))
end

figure(5)
%surf(pXi)
%surf(exp(pXi), 'EdgeColor', 'none', 'FaceColor', 'interp')
surf(resid, 'EdgeColor', 'none', 'FaceColor', 'interp')
colormap(jet)
view([0 0])
colorbar

figure(6)
i=20;
plot3(x(:,i),x(:,i+1),x(:,i+2))
xlabel('x1')
ylabel('x2')
zlabel('x3')
break

Theta = poolData(x,r,1,0);
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end


%A = Xi(2:r+1,1:r-1)';
A = pXi(2:r+1,1:r-1)';
B = A(:,r);
A = A(:,1:r-1);
%
L = 1:50000;
sys = ss(A,B,eye(r-1),0*B);

X0=x(1,1:r-1);         % initial state vector
Forc=x(L,r);           % Forcing term
[y,t] = lsim(sys,Forc,dt*(L-1),X0);

%% SAVE DATA (OPTIONAL)
% save ./DATA/FIG01_LORENZ.mat
