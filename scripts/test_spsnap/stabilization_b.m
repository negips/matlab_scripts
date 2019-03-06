%     Test Boostconv and spsnap

clear
clc
%close all

lafs = 20;

rng('default');

n = 50;    % Matrix Size

% No of eigenvalues greater than 1.
np = 20;
% How unstable do we want the eigenvalues?
uns_range = 0.01;

theta = 4*pi/9*rand(np,1);
e1  = (cos(theta) + uns_range*rand(np,1)) + 1i*sin(theta);
uns_e = [e1; conj(e1)];
% Set first eigenvector = 1 signifying beseflow

ne2 = n-length(e1);
st_range = 0.5;
theta = 4*pi/9*rand(ne2,1);
damp = 0.80 + 0.18*rand(ne2,1);
e2  = damp.*cos(theta) + 1i*sin(theta);
stb_e = [e2; conj(e2)];

Aeig = [1.0; uns_e; stb_e];
lr = real(Aeig);
li = imag(Aeig);

neig = length(Aeig);    % This is my matrix size
                        % Hopefully no repeated eigenvalues

figure(1);
scatter(lr,li); hold on
scatter(real(uns_e),imag(uns_e), 'r')
theta = linspace(0,2*pi,1000);
plot(cos(theta),sin(theta), '--k')
xlabel('$\lambda_{r}$', 'FontSize', lafs)
ylabel('$\lambda_{i}$', 'FontSize', lafs)


B = rand(neig);
B = B + B';

bnorm = norm(B);
disp(['Norm(B)=', num2str(bnorm)])
B = B/bnorm;
Binv = inv(B);

A = Binv*diag(Aeig)*B;

anorm = norm(A);
disp(['Norm(A)=', num2str(anorm)])

x0 = rand(neig,1) + 1i*rand(neig,1);

r = x0;
rnorm = norm(x0);

% This is the RHS
b = rand(neig,1) + 1i*rand(neig,1);

% SPSNAP parameters
ifsnap = 1;
skryl = 50;
sfreq = 1;
ifinit = 0;
X   = zeros(neig,skryl);
Y   = zeros(neig,skryl);
Z  = zeros(neig,skryl);
W  = zeros(neig,skryl);
vold  = zeros(neig,1);
vol1 = zeros(neig,1);
%

ik = 0;     % Subspace counter
vlen = neig;
resio = 1;
plotio  = 100;

niters = 2000;
residuals = zeros(niters,1);
x0norm    = zeros(niters,1);

for i=1:niters

   xi = A*x0;

   rnorm = norm(xi-x0);

   [xi_o,x0_o,X_o,Y_o,Z_o,W_o,vol1_o,vold_o,rnorm_o,ifinit_o,ik_o] = SpSnapOrtho_b(xi,x0,b,i,vol1,vold,X,Y,Z,W,ifsnap,sfreq,ifinit,ik,skryl,vlen);
   dv_o = [];
   ksize_o = 0;

   xi=xi_o;
   dv=dv_o;
   X=X_o;
   Y=Y_o;
   Z=Z_o;
   W=W_o;
   vol1=vol1_o;
   vold=vold_o;
   rnorm=rnorm_o;
   ifinit=ifinit_o;
   ik=ik_o;
   x0=x0_o;
   ksize = ksize_o;

   if mod(i,resio)==0
     disp([num2str(i),' Residual=',num2str(rnorm)])
     pause(0.001)
   end  

   residuals(i)=rnorm;
   x0norm(i) = norm(x0);

   if mod(i,plotio)==0
     figure(4)
     plot(real(x0))
     pause(0.01)
   end

   if i>10 && rnorm<1e-20
     break
   end  
  

end   

figure(2)
semilogy(residuals(1:i)); hold on
semilogy(x0norm(1:i)); hold on

figure(4)
plot(real(xi))


