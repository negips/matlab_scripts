%     Test Boostconv and spsnap

clear
clc
close all

lafs = 20;

%rng('default');

n = 50;    % Matrix Size

% No of eigenvalues greater than 1.
np = 1;
% How unstable do we want the eigenvalues?
uns_range = 0.0001;

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

bkryl = 15;
vin   = zeros(neig,bkryl);
vout  = zeros(neig,bkryl);
vold  = zeros(neig,1);
vol1 = zeros(neig,1);

ifboost = 1;
bfreq = 10;
boostio = 1;
plotio  = 100;
ib = 0;     % BoostConv counter
zro = zeros(neig,1);
c   = zeros(bkryl,1);

niters = 5000;
residuals = zeros(niters,1);
for i=1:niters

   xi = A*x0;

   rnorm = norm(xi-x0);
   x0 = xi;

   if (ifboost)
     if i==1
       vol1 = zro;
       vold  = zro;
        
       vold  = xi;
       ifinit = 0;
     end
   end


%  Apply Boostconv
   if (ifboost)
     if mod(i,bfreq)==0
       if (~ifinit)
          ib = 1;
          dv = (xi - vol1);

          vin(:,1)=dv;
          vout(:,1)=dv;
          ifinit = 1;
          
          vol1  = xi;        % not in code. Probably initialization bug
       else
          
          dv = (xi - vol1);

          vout(:,ib) = vout(:,ib)-dv;             % vout = rn_1 - rn
          vin(:,ib)  = vin(:,ib)-vout(:,ib);      % vin  = rn

          c_bc  = vout(:,1:ib)'*dv;
          dd_bc = vout(:,1:ib)'*vout(:,1:ib);
          c2    = dd_bc\c_bc;               % invert dd_bc
          c(1:length(c2),1) = c2;

          ib = mod(ib,bkryl)+1;
          vout(:,ib) = dv;                  % vout+1 = rn

%         New residual
          dv = dv + vin*c;
          vin(:,ib) = dv;

          rnorm = norm(dv);
         
          x0 = vol1+dv;
          vol1 = x0;
          vold = x0;            % also not in code

       end
     else
       vold = (xi-vold);
       rnorm = norm(vold);
       vold = xi;

     end          % if mod(i,bfreq)
   end            % ifboost     

   if mod(i,boostio)==0
     disp([num2str(i),' Residual=',num2str(rnorm)])
     pause(0.01)
   end  

   residuals(i)=rnorm;

   if mod(i,plotio)==0
     figure(3)
     plot(real(xi))
     pause(0.01)
   end  
  

end   

figure(2)
plot(residuals(1:i))


