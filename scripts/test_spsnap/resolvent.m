% Calculate Resolvent using Matrix Logarithm

clear
clc
close all

genmatrix

DT = 1.0;
[L,flag] = logm(A);
L=1/DT*L;

% ee = eig(L);
% lrL = real(ee);
% liL = imag(ee);
% 
% figure(2);
% scatter(lrL,liL, '.'); hold on
% xlabel('$\lambda_{r}$', 'FontSize', lafs)
% ylabel('$\lambda_{i}$', 'FontSize', lafs)
% title('Operator: $L$')
% 
% L2=expm(DT*L);
% ee = eig(L2);
% lrL2 = real(ee);
% liL2 = imag(ee);
% 
% figure(1)
% scatter(lrL2,liL2, '.k'); hold on
% title('Operator: $e^{tL}$')

% Assume random forcing function
Omega = 0.5;

rng('default')
fcap = rand(neig,1);
aa = linspace(0,2*pi,neig)';
%fcap = sin(aa);

II = eye(neig);
io = sqrt(-1);

ResolventOp = (io*Omega*II - L);
response  = ResolventOp\fcap; 

figure(3)
plot(response, 'o', 'LineWidth', 2); hold on
ylabel('$\hat{q}$')

% Matrixlog method

DT = 0.1;
Lnew = (L - io*Omega*II);

% New time-stepper
Anew = expm(DT*Lnew);

ee = eig(Anew);
lrAnew = real(ee);
liAnew = imag(ee);

figure(2);
scatter(lrAnew,liAnew, '.'); hold on
theta = linspace(0,2*pi,1000);
plot(cos(theta),sin(theta), '--k')
xlabel('$\lambda_{r}$', 'FontSize', lafs)
ylabel('$\lambda_{i}$', 'FontSize', lafs)
title('Operator: $e^{tL''}$')



rhs = DT*fcap;

b = rhs/norm(rhs);

exact_logAnew_u = logm(Anew)*b;


% Arnoldi to calculate log(Anew)*u
nk = neig;
Qk   = zeros(neig,nk+1);
Hess = zeros(nk+1,nk);
Qk(:,1) = b;
u = b;
resid = 1;
tol=1.0e-10;
resnorm = zeros(nk,1);
ngs = 2;
for k=1:nk+1
  u = Anew*u;
  h=zeros(k,1);
  for ig = 1:ngs
    g = Qk(:,1:k)'*u;
    u = u - Qk(:,1:k)*g;
    h = h + g;
  end 
  uortho = u;
%  uortho = u - Qk(:,1:k)*h;
  beta = norm(uortho);
  u = uortho/beta;
  Qk(:,k+1)=u;
  Hess(1:k+1,k) = [h; beta];
  H = Hess(1:k,1:k);
% Since norm of b=1
% Log(Anew)*u ~ Qk*Log(H)*e1
  e1=zeros(k,1);
  e1(1)=1;
  logAnew_u = Qk(:,1:k)*logm(H)*e1;
  resid = norm(exact_logAnew_u - logAnew_u);
  resnorm(k)=resid;    

  if (k>1)
    disp([num2str(k) ' Log(A)u Residual norm: ', num2str(resid,5)])
    logAnew_old = logAnew_u;
  else
    logAnew_old = logAnew_u;
  end  
  if resid<tol
    break
  end  
end   
 
figure(4)
plot(exact_logAnew_u, 'o', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(logAnew_u, '.', 'LineWidth', 2, 'MarkerSize', 12); hold on
ylabel('$log(Anew)b$')


figure(5)
semilogy(resnorm(1:k))














