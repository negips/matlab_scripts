% SNAP as implemented in paper

clear
clc
close all

genmatrix

x0 = rand(neig,1) + 1i*rand(neig,1);

r = x0;
rnorm = norm(x0);

% This is the RHS
b = rand(neig,1) + 1i*rand(neig,1);

% SPSNAP parameters
ifsnap = 1;
skryl =  25;
sfreq = 1;
ifinit = 0;
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

ph = 1;
mode = 0;
bb = b'*b;
p1gmres = 10;
p2gmres = 5;

Kryl_x = zeros(neig,p1gmres);
Kryl_p1=zeros(neig,p1gmres);

Kryl_x2 = zeros(neig,p2gmres);
Kryl_p2 = zeros(neig,p2gmres);

nnull = 50;
X = zeros(neig,nnull);
AX  = zeros(neig,nnull);
AbX = zeros(neig,nnull);

for i=1:niters

   xi = A*x0;

   rnorm = norm(xi-x0);

   if (ph==1)   
     if mode==0               % Skip first vector
      x0 = xi;
      mode=1;
      continue
     elseif mode==1
       v0 = xi;
       x0 = v0;
       mode=2;
       continue
     elseif mode==2
       v1 = -(xi - b*(b'*xi)/bb);
       x0 = v1;
       mode=3;
       Kryl_x(:,1) = x0;
       ik=0;
       continue
     elseif mode==3
       ik=ik+1;
       t1 = xi - b*(b'*xi)/bb;
       Kryl_p1(:,ik)=t1;
       if ik==p1gmres
         ktk = Kryl_p1'*Kryl_p1;
         ktv = -Kryl_p1'*v1;
         t_red = ktk\ktv;
         t = Kryl_x*t_red;
         r = Kryl_p1*t_red - v1;
         w1 = v0+t;           % First approximate null space
         X(:,1) = w1;
         inull = 1;
         x0 = w1;
         ph = 2;              % go to phase 2
         mode=0;
         continue
       else 
         Kryl_x(:,ik) = xi;
         x0           = xi;
         continue
       end
     else
       disp(['Unknwn mode', num2str(mode)])
       break
     end
     
   else           % Phase 2
     if mode==0
       wsq = w1'*w1;
       epsi = xi - b*(b'*xi)/bb;

       AX(:,1)  = xi;
       AbX(:,1) = epsi;
       v1 = -(epsi - w1*(w1'*epsi)/wsq);
       x0 = v1;
       Kryl_x2(:,1) = x0;
       ik = 0;
       mode=1;
       continue
     elseif mode==1
       ik=ik+1;
       t1 = xi - w1*(w1'*xi)/wsq;
       Kryl_p2(:,ik)=t1;
       if ik==p2gmres
         ktk = Kryl_p2'*Kryl_p2;
         ktv = -Kryl_p2'*v1;
         t_red = ktk\ktv;
         t = Kryl_x2*t_red;
         r = Kryl_p2*t_red + v1;
         w2 = w1+t;           % Next approximate null space vector
         inull = inull+1;
         X(:,inull) = w2;
         x0 = w2;
         mode=3;
         continue
       else 
         Kryl_x2(:,ik) = xi;
         x0            = xi;
         continue
       end
     elseif mode==3  
       t = xi - b*(b'*xi)/bb;
       AX(:,inull) = xi;
       AbX(:,inull) = t;
       [U,S,V] = svd(AbX(:,1:inull),'econ');
       s=diag(S);
       smin = s(end);
       vmin = V(:,end);
       w1   = X(:,1:inull)*vmin;
       epsi = smin*U(:,end);

       beta = bb/(b'*AX(:,1:inull)*vmin);
       rnorm = abs(beta)*norm(epsi);

       soln = beta*w1;

%       disp(['Sigma0=', num2str(smin)])
       disp(['beta=', num2str(abs(beta))])
%       disp(['Residual=', num2str(rnorm)])
     

       wsq = w1'*w1;
       v1  = -(epsi - w1*(w1'*epsi)/wsq);
       x0  = v1;
       Kryl_x2(:,1) = x0;
       ik = 0;
       mode=1;
       if inull==nnull
         break
       end 
     end  

   end  


%   if mod(i,resio)==0
%     disp([num2str(i),' Residual=',num2str(rnorm)])
%     pause(0.001)
%   end  

%   residuals(i)=rnorm;
%   x0norm(i) = norm(x0);
  

end   









