% This script runs the temporal 1D stability code
%
% (c) Ardeshir Hanifi 2014
%
clear all
clc
close all
%% Numerical parameters

N = 200; % number of Cehebychev modes

%% Select the collocation nodes and set up derivative operator
M=4; % Highest order of derivative required
[y, D] = chebdif(N, M);

%% Velocity profile corresponding to Backflow (Alam and Sandham 2000)
ymax = 2;
Jac=ymax/2;
%y1 = (y+1)*ymax/2;
y1 = y;
a1=1.16;
b1=1.55;
U.u  = tanh(y1) - 2*a1*tanh(y1/b1)./(cosh(y1/b1).^2);

dstar = trapz(y1,(1 - U.u));
dstar = abs(dstar);
ynorm = 1;

y1 = y1/ynorm;
ymax = max(y1);
figure(100)

D(:,:,1) = D(:,:,1)/(Jac);
D(:,:,2) = D(:,:,2)/(Jac^2);
D(:,:,3) = D(:,:,3)/(Jac^3);
D(:,:,4) = D(:,:,4)/(Jac^4);
D1=D(:,:,1);
D2=D(:,:,2);

U.uy = D1*U.u;
U.uyy= D2*U.u;
Re=500*dstar;                 % Redstar

% Plane Poiseuille flow
U.u=1-y.^2;
U.uy=-2*y;
U.uyy=-2+y*0;
Re=2000;


% Testing with Couette
%U.u = y;
%U.uy = D1*U.u;
%U.uyy = D2*U.u;

plot(y1,U.u); hold on
plot(y1,U.uy, 'r')
plot(y1,U.uyy, 'k')
grid on

%% Flow parameters
beta=0.00;

omega_min=0.0
omega_max=0.0;
nomega=1;
omega_range = linspace(omega_min,omega_max,nomega);

cols = lines(nomega);

%% Set up the eigenvalue problem

%uns_ar = [];
%uns_ai = [];
%uns_wr = [];
%uns_wi = [];
%uns_ind = [];
nuns = 0;

nfamilies=2;
ifverify=1;
cols_fam = lines(nfamilies);

neigs = 0;
zi = sqrt(-1);

omega_r = 0.3;
for ios=1:nomega
  omega = (omega_r + zi*omega_range(ios));
  [A,B] = OS_spatial_extended(U, D, omega, beta, Re);
  
  %% Solve the eigenvalue problem
  [eigvec, eigval] = eig(A,B);
  eigval=diag(eigval)*(-1);
  
  %% Plot spectrum and eigenvectors
%  header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
%                {', beta= '},{num2str(beta)});
  header=strcat({'Re= '},{num2str(Re)},{', beta= '},{num2str(beta)});

  %plot_OS(eigval,eigvec,y1,header)

  ind1 = find(abs(imag(eigval))<1000000);
  lambda =  eigval(ind1);
  lr = real(lambda);
  li = imag(lambda);
  evec2 = eigvec(N+1:2*N,ind1); 
 
  figure(1)
  plot(lr,li,'.', 'Color', cols(ios,:)); hold on
  xlims=get(gca,'XLim');
  grid on
  xlim([-0.5 1.5])
  ylim([-0.1 1]) 
%  plot(xlims,[0 0], '--')
%  xlim([-20,20]);
  title(header);
  ylabel('$imag(\alpha)$')
  xlabel('$real(\alpha)$');

   if neigs==0
%     [val ind] = max(imag(eigval));

     for ii=1:nfamilies 
       found=0; 
       while ~found 
         figure(1)   
         [xp,yp,button] = ginput(1);
         a=xp+sqrt(-1)*yp;
         [c,ind]=min(abs(lambda-a));
         b=lambda(ind);   
         if ifverify 
           figure(10)
           plot(y1,real(evec2(:,ind).*exp(-b*y)), 'b','Marker', '.'); hold on
           plot(y1,imag(evec2(:,ind).*exp(-b*y)), 'r', 'Marker', '.'); hold on
           grid on 

           userin = input('Select Mode (1/0)?');
           if userin
              found=1;
           end
           close(10)   
         else
           found=1;
         end      % verify
       end        % ~found
       % Mark selected
       figure(1)
       plot(lambda(ind), 'ok', 'MarkerSize', 8)      

       ref_lambda{ii} = lambda(ind);

       uns_ar{ii} = real(lambda(ind));
       uns_ai{ii} = imag(lambda(ind));
       uns_wr{ii} = real(omega);
       uns_wi{ii} = imag(omega);
       uns_ind{ii} = ind;
       eigold{ii} = lambda(ind);
     end          % nfamilies       
     neigs=neigs+1; 
   else

     for ii=1:nfamilies
       [c,ind]=min(abs(lambda-eigold{ii}));
       uns_ar{ii} = [uns_ar{ii} real(lambda(ind))];
       uns_ai{ii} = [uns_ai{ii} imag(lambda(ind))];
       uns_wr{ii} = [uns_wr{ii} real(omega)];
       uns_wi{ii} = [uns_wi{ii} imag(omega)];
       uns_ind{ii} = [uns_ind{ii} ind];
       eigold{ii} = lambda(ind);
     end 
     neigs=neigs+1;

   end      % neigs==0

   nuns=nuns+1;
%   figure(2)
%   plot(alpha,real(eigval(ind)), '.', 'Color', cols(ios,:), 'MarkerSize', 9); hold on

end

if (nuns>1)
  figure(2)
  for ii=1:nfamilies    
    plot(uns_ar{ii},uns_ai{ii}, 'Color', cols_fam(ii,:)); hold on
    grid on
  end    
%
%  cg = gradient(uns_wr,uns_ar);
%  figure(3)
%  plot(uns_ar,cg)
%  grid on     
end
