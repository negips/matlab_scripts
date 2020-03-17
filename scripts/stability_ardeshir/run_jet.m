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

%% Velocity profile corresponding to Plane Poiseuille Flow
ymax = 5;
y1 = (y)*ymax;
U.u  = sech(y1).^2;


dstar = trapz(y1,(1 - U.u));
dstar = abs(dstar);
ynorm = 1;

y1 = y1/ynorm;
ymax = max(y1);
%figure(10)
%plot(y1,U.u)


D(:,:,1) = D(:,:,1)/(ymax);
D(:,:,2) = D(:,:,2)/(ymax^2);
D(:,:,3) = D(:,:,3)/(ymax^3);
D(:,:,4) = D(:,:,4)/(ymax^4);

D1=D(:,:,1);
D2=D(:,:,2);


U.uy = D1*U.u;
U.uyy= D2*U.u;

Re=10;                 % Redstar
%% Flow parameters
alpha=0.5;
beta=0.000;

alpha_min=0.1;
alpha_max=3.5;
nalpha=50;
alpha_range = linspace(alpha_min,alpha_max,nalpha);
%nalpha=1;
%alpha_range=[0.45];

cols = lines(nalpha);

%% Set up the eigenvalue problem

uns_wr = [];
uns_wi = [];
uns_ar = [];
uns_ind = [];
nuns = 0;

neigs = 0;
for ios=1:nalpha
  alpha=alpha_range(ios);
  [A,B] = OS_temp(U, D, alpha, beta, Re);
  
  %% Solve the eigenvalue problem
  [eigvec, eigval] = eig(A,B);
  eigval=diag(eigval);
  
  %% Plot spectrum and eigenvectors
%  header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
%                {', beta= '},{num2str(beta)});
  header=strcat({'Re= '},{num2str(Re)},{', beta= '},{num2str(beta)});

  %plot_OS(eigval,eigvec,y1,header)
  
  figure(1)
  plot(real(eigval),imag(eigval),'.', 'Color', cols(ios,:)); hold on
  ylim([-2,2]);
  xlims=get(gca,'XLim');
  plot(xlims,[0 0], '--')
  title(header);
  ylabel('$\omega_{i}$','interpreter','latex')
  xlabel('$\omega_{r}$','interpreter','latex')

  [val ind] = max(imag(eigval));
%  if (val>0)
%    nuns = nuns + 1;  
%    uns_wr = [uns_wr real(eigval(ind))];
%    uns_wi = [uns_wi imag(eigval(ind))];
%    uns_ar = [uns_ar alpha];
%    uns_ind = [uns_ind ind];  
%  
%    figure(2)
%    plot(alpha,real(eigval(ind)), '.', 'Color', cols(ios,:), 'MarkerSize', 9); hold on
%  end

   if neigs==0
%     [val ind] = max(imag(eigval)); 
     [xp,yp,button] = ginput(1);
     a=xp+sqrt(-1)*yp;
     [c,ind]=min(abs(eigval-a));
     disp(['(',num2str(real(eigval(ind))), ',' , num2str(imag(eigval(ind))), 'i)'])
     uns_wr = [uns_wr real(eigval(ind))];
     uns_wi = [uns_wi imag(eigval(ind))];
     uns_ar = [uns_ar alpha];
     uns_ind = [uns_ind ind];
     eigold = eigval(ind);
     neigs=neigs+1; 
   else
     [c,ind]=min(abs(eigval-eigold));
     uns_wr = [uns_wr real(eigval(ind))];
     uns_wi = [uns_wi imag(eigval(ind))];
     uns_ar = [uns_ar alpha];
     uns_ind = [uns_ind ind];
     eigold = eigval(ind);
     neigs=neigs+1; 
   end

   nuns=nuns+1;
   figure(2)
   plot(alpha,imag(eigval(ind)), '.', 'Color', cols(ios,:), 'MarkerSize', 9); hold on

end

if (nuns>1)
  figure(2)
  plot(uns_ar,uns_wi)
  xlabel('$\alpha_{r}$','interpreter','Latex')
  ylabel('$\omega_{i}$','interpreter','Latex')
  grid on

%  cg = gradient(uns_wr,uns_ar);
%  figure(3)
%  plot(uns_ar,cg)
%  grid on

  figure(1)
  plot(uns_wr,uns_wi)
  xlabel('$\omega_{r}$','interpreter','Latex')
  ylabel('$\omega_{i}$','interpreter','Latex')
  grid on
 
end
