function [dt dt2 dc dc2 lp cbc lpbc forc abcde D x wts] = MEFem(N,xs,xe,ifplot)

addpath '../';

jac = (xe-xs)/2;              % jacobian
dxbydxi = (xe-xs)/2;
dxibydx = 2/(xe-xs);

[x wts p]= lglnodes(N);
x =x(end:-1:1);
wts =wts(end:-1:1);

% build A
A1 = [];
for i = 0:N
     A1= [A1 x.^i];
end
abcde = zeros(N+1);

for i = 0:N
     % build b
     b = zeros(N+1,1);
     b(i+1) = 1;

     solns = A1\b;
     abcde(:,i+1) = solns;
end

D = [];
for i = 0:N
     j = (i-1);
     t1 = i*abcde(i+1,:);
     nans = isnan(t1);
     if max(nans)
          ind = find(nans);
          t1(ind) = 0;
     end
     D= [D; t1];
end

%abcde
%D

%% testing
xtemp = transpose(linspace(-1,1,200));
l2 = length(xtemp);
basis = zeros(l2,N+1);
for i = 0:N
val = zeros(l2,1);
for j = 0:N
     val = val + abcde(j+1,i+1)*xtemp.^j;
end
basis(:,i+1) = val;
end

if ifplot
  h1 = figure;
  hold on
  plot(xtemp,basis)
end

% Derivative function
deriv = zeros(l2,N+1);
for i = 0:N
val = zeros(l2,1);
for j = 0:N
  val1 = D(j+1,i+1)*(xtemp.^(j-1));
  if (j-1)<0
       val1 = zeros(l2,1);
  end
  val = val + val1;
end
deriv(:,i+1) = val;
end

if ifplot
  h2 = figure;
  hold on
  plot(xtemp,deriv)
end

% integral vector
dx_integral = zeros(2*N+1,1);
for j = 0:2*N
     xorder = j;
     dx_integral(j+1) = (1^(xorder+1) - (-1)^(xorder+1))/(xorder+1);
end

% matrix in front of d/dt
dt = zeros(N+1);

for j = 0:N

Lj = abcde(:,j+1);

xj = x(j+1);
x_p_vec=zeros(N+1,1);
for i = 0:N
    x_p_vec(i+1) = xj^i;
end
Lj_xj = Lj.*x_p_vec; 
% [x_p_vec Lj Lj_xj]

integral = jac*wts(j+1)*sum(Lj_xj)*sum(Lj_xj);
dt(j+1,j+1) = integral;

end

%dt

%% Analytically compute integral
dt2 = zeros(N+1);

for j = 0:N

Lj = abcde(:,j+1);

integral = 0;
xj = x(j+1);
for i = 0:N
for k = 0:N
     xorder=k+i;
     coeff = Lj(i+1)*Lj(k+1);
     integral = integral + coeff*dx_integral(xorder+1);
end
end

dt2(j+1,j+1) = jac*integral;

end


% matrix for vdu/dx           % variational form: -udv/dx

dc = zeros(N+1);

for j = 0:N      % Test function (v)

     DLj = D(:,j+1);

for i = 0:N      % Trial function (u)

     Li = abcde(:,i+1);
     x_pi_vec=zeros(N+1,1);
     xi = x(i+1);
     for k = 0:N                   
         x_pi_vec(k+1) = xi^k;
     end
     x_pj_vec=zeros(N+1,1);
     for k = 1:N                   
         x_pj_vec(k+1) = xi^(k-1);
     end

     Li_xj = Li.*x_pi_vec;
     dLj_xj = DLj.*x_pj_vec; 

%     if j == 2
%          [x_p_vec dLi dLi_xj]
%     end
     integral = jac*dxibydx*wts(i+1)*sum(Li_xj)*sum(dLj_xj);

     dc(j+1,i+1) = integral;

end

end

dc = -dc;

%% Boundary term  for convection        % after integration by parts: uv(1) - uv(-1)
cbc = zeros(N+1);
loop = 1;

%

dc2 = zeros(N+1);                       % convection operator without integration by parts:    vdu/dx

for j = 0:N      % Test function (v)

     Lj = abcde(:,j+1);

     xj = x(j+1);
     x_p_vec=zeros(N+1,1);
     for k = 0:N
         x_p_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_p_vec; 
     % [x_p_vec Lj Lj_xj]

     integral = 0;
for i = 0:N      % Trial function (u)

     dLi = D(:,i+1);
     x_p_vec=zeros(N+1,1);
     for k = 1:N                   % k=0 ==> has coefficient of 0. Hence skipped.
         x_p_vec(k+1) = xj^(k-1);
     end
     dLi_xj = dLi.*x_p_vec;
%     if j == 2
%          [x_p_vec dLi dLi_xj]
%     end
     integral = jac*dxibydx*wts(j+1)*sum(Lj_xj)*sum(dLi_xj);

     dc2(j+1,i+1) = integral;

end

end

dc2 = dc2;

%% Laplacian term                  % Integration by parts:  dudx.dvdx
lp = zeros(N+1);

for j = 0:N                   % Test function (v)
     dLj = D(:,j+1);
for i = 0:N                   % Trial function (u)
     dLi = D(:,i+1);
     
     integral = 0;
for k = 0:N

     xk = x(k+1);

     x_p_vec=zeros(N+1,1);
     for kk = 1:N                   % k=0 ==> has coefficient of 0. Hence skipped.
         x_p_vec(kk+1) = xk^(kk-1);
     end

     dLi_xk = dLi.*x_p_vec;
     dLj_xk = dLj.*x_p_vec; 
     
     integral = integral + (jac*(dxibydx)^2)*wts(k+1)*sum(dLj_xk)*sum(dLi_xk);

end

     lp(j+1,i+1) = integral;

end

end

lp = -lp;


lpbc = zeros(N+1);             % Boundary terms for laplacian operator. vdudx(1) - vdudx(-1)
loop = 1;
for j = 0:loop      % Test function (v)
     
     Lj = abcde(:,j*N+1);

     xj = x(j*N+1);
     x_pj_vec=zeros(N+1,1);
     for k = 0:N
         x_pj_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_pj_vec; 
     % [x_p_vec Lj Lj_xj]

for i = 0:N      % Trial function (u)

     dLi = D(:,i+1);
     x_pi_vec=zeros(N+1,1);
     for k = 1:N                   % derivative. hence k=0 skipped.
         x_pi_vec(k+1) = xj^(k-1);
     end
     dLi_xj = dLi.*x_pi_vec;

     integral = jac*dxibydx*((-1)^j)*sum(dLi_xj)*sum(Lj_xj);

     lpbc(j*N+1,i+1) = integral;

end

end

lpbc = lpbc;



% Hack for testing ...
%for i = 1:N+1
%     dc2(i,i) = 0;
%end

%%

%% Build forcing

%% Build basis change matrix
%% Taken directly from nek
nx = N+1;
kj = 0;
n = nx-1;
for j = 1:nx
     z=x(j);
     Lj = legendrePoly(n,z);
     kj = kj+1;
     pht(kj) = Lj(1);
     kj = kj+1;
     pht(kj) = Lj(2);
     for k = 3:nx
          kj = kj+1;
          pht(kj) = Lj(k);% - Lj(k-2);
     end
end

pht = reshape(pht,nx,nx);
spectonodal = transpose(pht);
phi = pht;

boyd = 0;
if boyd
  for i = 3:nx
    for j = 1:nx
      phi(i,j) = pht(i,j) - pht(i-2,j);
    end
  end
end
boydstonodal = transpose(phi);

%%
% matrix for filter forcing

G = zeros(N+1);
G(N+1,N+1)=1;
%G(N,N) = 0.5;
%G(N,N) =1;
%G(N-1,N-1) =0.1;

%G = eye(N+1);

fil_mat = boydstonodal*G*inv(boydstonodal);

fil = zeros(N+1);

for j = 0:N      % Test function (v)
  Lj = abcde(:,j+1);

  xj = x(j+1);
  x_p_vec=zeros(N+1,1);
  for k = 0:N
      x_p_vec(k+1) = xj^k;
  end
  Lj_xj = Lj.*x_p_vec; 

  Fi = fil_mat(j+1,:);
  
  integral = jac*wts(j+1)*sum(Lj_xj).*Fi;

  forc(j+1,:) = integral;
end


end
