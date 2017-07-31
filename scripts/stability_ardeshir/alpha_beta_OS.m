% This script runs the temporal 1D stability code
%
% (c) Ardeshir Hanifi 2014
%
clear all
clc
close all
%% Numerical parameters

N = 151; % number of Cehebychev modes

%% Select the collocation nodes and set up derivative operator
M=4; % Highest order of derivative required
[y, D] = chebdif(N, M);

%% Velocity profile corresponding to Plane Poiseuille Flow
ymax = 16;
y1 = (y+1)*ymax/2;
a1=1.16;
b1=1.55;
%U.u=1-y.^2;
%U.uy=-2*y;
%U.uyy=-2+y*0;
U.u  = tanh(y1) - 2*a1*tanh(y1/b1)./(cosh(y1/b1).^2);

dstar = trapz(y1,(1 - U.u));
dstar = abs(dstar);
ynorm = 1;

y1 = y1/ynorm;
ymax = max(y1);
figure(10)
plot(y1,U.u)


D(:,:,1) = D(:,:,1)/(ymax)*2;
D(:,:,2) = D(:,:,2)/(ymax^2)*(2^2);
D(:,:,3) = D(:,:,3)/(ymax^3)*(2^3);
D(:,:,4) = D(:,:,4)/(ymax^4)*(2^4);

D1=D(:,:,1);
D2=D(:,:,2);


U.uy = D1*U.u;
U.uyy= D2*U.u;

Re=500;                 % Redstar
%% Flow parameters
alpha=1.5;
beta=0.000;

alpha_min=-0.151;
alpha_max=0.152;
nalpha=150;
alpha_range = linspace(alpha_min,alpha_max,nalpha);

beta_min=-0.1;
beta_max=0.1;
nbeta=100;
beta_range = linspace(beta_min,beta_max,nbeta);

cols = lines(nalpha*nbeta);

%% Set up the eigenvalue problem

uns_wr = zeros(nbeta,nalpha);
uns_wi = zeros(nbeta,nalpha);
uns_ar = zeros(nbeta,nalpha);
uns_br = zeros(nbeta,nalpha);

nuns = 0;
neigs = 0;
nplots = 0;
ifplot = 1;

for ios=1:nalpha
  for jos=1:nbeta

    if (nplots>0)
%      pause(0.001)
%      delete(eigplot)
    end    

    alpha=alpha_range(ios);
    beta=beta_range(jos);  
    [A,B] = OS_temp(U, D, alpha, beta, Re);
    
    %% Solve the eigenvalue problem
    [eigvec, eigval] = eig(A,B);
    eigval=diag(eigval);
    
    %% Plot spectrum and eigenvectors
    header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
                  {', beta= '},{num2str(beta)});
%    header=strcat({'Re= '},{num2str(Re)},{', beta= '},{num2str(beta)});

    %plot_OS(eigval,eigvec,y1,header)

    if nplots==0
      
      figure(1)
      eigplot = plot(real(eigval),imag(eigval),'.', 'Color', cols(ios,:)); hold on
      xlims=get(gca,'XLim');
%      plot(xlims,[0 0], '--')
      ylim([-1,0.2]);
      title(header);
      ylabel('imag(\omega)')
      xlabel('real(\omega)')
      nplots=nplots+1;
    end    

    if neigs==0
%      [val ind] = max(imag(eigval)); 
      [xp,yp,button] = ginput(1);
      a=xp+sqrt(-1)*yp;
      [c,ind]=min(abs(eigval-a));
      uns_wr(jos,ios) = real(eigval(ind));
      uns_wi(jos,ios) = imag(eigval(ind));
      uns_ar(jos,ios) = alpha;
      uns_br(jos,ios) = beta;

      eigold = eigval(ind);

      neigs=neigs+1; 
    else
      [c,ind]=min(abs(eigval-eigold));
      uns_wr(jos,ios) = real(eigval(ind));
      uns_wi(jos,ios) = imag(eigval(ind));
      uns_ar(jos,ios) = alpha;
      uns_br(jos,ios) = beta;

      if (jos==nbeta) && (nbeta>1)
        eigold = uns_wr(1,ios) + sqrt(-1)*uns_wi(1,ios);
      else
        eigold = eigval(ind);
      end

      neigs=neigs+1; 
    end

  end     

  if ifplot
    figure(2)
    branchplot = plot(uns_wr(:,ios),uns_wi(:,ios),'.', 'Color', cols(ios,:)); hold on
    header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)}); 
    title(header);
    ylabel('imag(\omega)')
    xlabel('real(\omega)')
  end  

end

if (neigs>1)
  figure(3)
  surf(uns_ar,uns_br,uns_wr, 'EdgeColor', 'None', 'FaceColor', 'interp')
  xlabel('\alpha')
  ylabel('\beta')
  zlabel('\omega_{r}')             

  [Cgx Cgy] = gradient(uns_wr,alpha_range,beta_range);
  Cg2 = Cgx.*Cgy;   
  figure(4)
  surf(uns_ar,uns_br,Cg2, 'EdgeColor', 'None', 'FaceColor', 'interp')
  xlabel('\alpha')
  ylabel('\beta')
  zlabel('d\omega/d\alpha')             

end


