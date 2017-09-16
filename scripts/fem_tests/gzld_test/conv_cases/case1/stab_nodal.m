% testing stuff...

%close all
%clear
%clc

% Taking a periodic domain
% ==> U0 = Un;

%addpath '../';

N = 5;

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

abcde;
D;

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

h1 = figure;
hold on
plot(xtemp,basis)

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

h2 = figure;
hold on
plot(xtemp,deriv)


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

integral = wts(j+1)*sum(Lj_xj)*sum(Lj_xj);
dt(j+1,j+1) = integral;

end

dt;

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

dt2(j+1,j+1) = integral;


end

%dt2

%% small test
rat = zeros(N+1);
for i = 0:N
     rat(i+1,i+1) = dt(i+1,i+1)/dt2(i+1,i+1);
end
%rat

% matrix for cvdu/dx
conv = 1.0;
dc = zeros(N+1);

for j = 0:N      % Test function (v)

     DLj = D(:,j+1);

for i = 0:N      % Trial function (u)

     Li = abcde(:,i+1);
     x_pi_vec=zeros(N+1,1);
     x_pj_vec=zeros(N+1,1);
     xi = x(i+1);
     for k = 0:N                   
         x_pi_vec(k+1) = xi^k;
     end
     for k = 1:N                   
         x_pj_vec(k+1) = xi^(k-1);
     end


     Li_xj = Li.*x_pi_vec;
     dLj_xj = DLj.*x_pj_vec; 

%     if j == 2
%          [x_p_vec dLi dLi_xj]
%     end
     integral = wts(i+1)*sum(Li_xj)*sum(dLj_xj);

%     if i == N+1
%          dc(j+1,1) = dc(j+1,1) + integral;
%     else
          dc(j+1,i+1) = integral;
%     end

end

end

dc = -conv*dc;

%% Boundary term
bc = zeros(N+1);

for j = 0:1      % Test function (v)
     
     Lj = ((-1)^j)*abcde(:,j*N+1);

     xj = x(j*N+1);
     x_pj_vec=zeros(N+1,1);
     for k = 0:N
         x_pj_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_pj_vec; 
     % [x_p_vec Lj Lj_xj]

for i = 0:N      % Trial function (u)

     Li = abcde(:,i+1);
     x_pi_vec=zeros(N+1,1);
     for k = 0:N                % k=0 ==> has coefficient of 0. Hence skipped.
         x_pi_vec(k+1) = xj^k;
     end
     Li_xj = Li.*x_pi_vec;

%     if j == 2
%          [x_p_vec dLi dLi_xj]
%     end
     integral = sum(Li_xj)*sum(Lj_xj);

%     if i == N+1
%          bc(j*N+1,1) = bc(j+1,1) + integral;
%     else
          bc(j*N+1,i+1) = integral;
%     end

end

end

bc = conv*bc;

dc2 = dc + bc;

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

boyd = 1;
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

%G = zeros(N+1);
%G(N+1,N+1) =1;
%G(N,N) = 0.5;
%G(N,N) =1;
%G(N-1,N-1) =1;

G = eye(N+1);

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
     
     integral = wts(j+1)*sum(Lj_xj).*Fi;

     forc(j+1,:) = integral;

end

chi = 0;
df_mat = chi*forc;

% du/dt = Au;
A_mat1 = (df_mat - dc2);
A_mat = inv(dt)*A_mat1;
e1 = eig(A_mat);

cmap = jet(2);

h3 = figure;
plot(real(e1),imag(e1), '*', 'MarkerSize', 12, 'Color', cmap(1,:));
grid on
title(num2str(chi), 'FontSize', 16);
hold on

A_mat = inv(dt2)*A_mat1;
e2 = eig(A_mat);
plot(real(e2),imag(e2), '*', 'MarkerSize', 12, 'Color', 'r');



