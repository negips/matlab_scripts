% testing stuff...

close all
clear
clc

N = 16;

[x wts p]= lglnodes(N);
x =x(end:-1:1)
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

abcde
D

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
dxi_integral = zeros(2*N+1,1);
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

dt
% matrix for cvdu/dx
conv = 1.0;
dc = zeros(N+1);

for j = 0:N      % Test function (v)

     Lj = abcde(:,j+1);

     xj = x(j+1);
     x_p_vec=zeros(N+1,1);
     for k = 0:N
         x_p_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_p_vec; 
     % [x_p_vec Lj Lj_xj]

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
     integral = wts(j+1)*sum(Lj_xj)*sum(dLi_xj);

     dc(j+1,i+1) = integral;

end

end

dc = conv*dc

%% Build forcing

deltat = 0.0002;
alpha_final = 1/deltat;
steps = 100;
alphas = -linspace(0,alpha_final,steps);

h3=figure;
for k = 1:steps 

clc
k

chi = alphas(k);
%alpha = -0;
G = zeros(N+1);
G(N+1,N+1) =1;
%G(N,N) = 0.5;
%G(N-1,N+1) =0.1;
%G = eye(N+1);
df_mat = chi*G*dt;

% du/dt = Au;
A_mat1 = (df_mat - dc);
A_mat = inv(dt)*A_mat1;
e = eig(A_mat);

lambdadt = e*deltat;

figure(h3)
plot(real(lambdadt),imag(lambdadt), '*', 'MarkerSize', 6);
grid on
title(num2str(chi), 'FontSize', 16);
hold on

pause(0.05)
end

