% This script runs the temporal 1D stability code
%
% (c) Ardeshir Hanifi 2014
%
clear all
clc
close all
%% Numerical parameters

N = 201; % number of Cehebychev modes

%% Select the collocation nodes and set up derivative operator
M=4; % Highest order of derivative required
[y, D] = chebdif(N, M);

%% Velocity profile corresponding to Plane Poiseuille Flow
ymax = 2;
y1 = y*ymax/2;
a1=1.16;
b1=1.55;
U.u=1-y.^2;
U.uy=-2*y;
U.uyy=-2+y*0;
%U.u  = tanh(y1) - 2*a1*tanh(y1/b1)./(cosh(y1/b1).^2);

dstar = trapz(y1,(1 - U.u));
dstar = abs(dstar);
ynorm = 1;

%y1 = y1/ynorm;
%ymax = max(y1);
%figure(10)
%plot(y1,U.u)

D(:,:,1) = D(:,:,1)/(ymax)*2;
D(:,:,2) = D(:,:,2)/(ymax^2)*(2^2);
D(:,:,3) = D(:,:,3)/(ymax^3)*(2^3);
D(:,:,4) = D(:,:,4)/(ymax^4)*(2^4);

D1=D(:,:,1);
D2=D(:,:,2);


U.uy = D1*U.u;
U.uyy= D2*U.u;

%Re=500*dstar;                 % Redstar
Re=1000;
%% Flow parameters
alpha=1.5;
beta=0.000;

alpha_min=1e-5;
alpha_max=1.2;
nalpha=200;
alpha_range = linspace(alpha_min,alpha_max,nalpha);

alpha_range = [1.0];
nalpha = 1;

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
  plot(real(eigval),imag(eigval),'.', 'Color', cols(ios,:), 'MarkerSize', 10); hold on
  xlims=get(gca,'XLim');
%  plot(xlims,[0 0], '--')
%  ylim([-1,0.2]);
  title(header);
  ylabel('$\omega_{i}$')
  xlabel('$\omega_{r}$')
  xlim([0 1])
  ylim([-1 0])

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

%   if neigs==0
%%     [val ind] = max(imag(eigval)); 
%     [xp,yp,button] = ginput(1);
%     a=xp+sqrt(-1)*yp;
%     [c,ind]=min(abs(eigval-a));
%     uns_wr = [uns_wr real(eigval(ind))];
%     uns_wi = [uns_wi imag(eigval(ind))];
%     uns_ar = [uns_ar alpha];
%     uns_ind = [uns_ind ind];
%     eigold = eigval(ind);
%     neigs=neigs+1; 
%   else
%     [c,ind]=min(abs(eigval-eigold));
%     uns_wr = [uns_wr real(eigval(ind))];
%     uns_wi = [uns_wi imag(eigval(ind))];
%     uns_ar = [uns_ar alpha];
%     uns_ind = [uns_ind ind];
%     eigold = eigval(ind);
%     neigs=neigs+1; 
%   end

   nuns=nuns+1;
   figure(2)
   plot(alpha,real(eigval(ind)), '.', 'Color', cols(ios,:), 'MarkerSize', 9); hold on

end

if (nuns>1)
  figure(2)
  plot(uns_ar,uns_wr)
  grid on

  cg = gradient(uns_wr,uns_ar);
  figure(3)
  plot(uns_ar,cg)
  grid on     
end
