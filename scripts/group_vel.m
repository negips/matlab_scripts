%  Generating a 2x2 system to figure out group velocity


% du/dt + a1*du/dx + a2*dv/dx -nu1*d2u/dx2 + a3*u = 0
% dv/dt + a4*du/dx + a5*dv/dx -nu2*d2v/dx2 + a6*u = 0

clear
clc
close all

% Using ansatz u=exp(kx-wt); v=exp(kx-wt);

k=1.5;
kall = -50:0.1:50;
nk = length(kall);

lr = zeros(nk,2);
li = zeros(nk,2);

for i=1:nk
  k  = kall(i);
  A  = BuildA(k);
  eg = eig(A);
  lr(i,:) = real(eg)';
  li(i,:) = imag(eg)';

end

plot(lr,li, '.', 'MarkerSize', 8)



%----------------------------------------------------------------------  
function A = BuildA(k)

      a1 = 1;
      a2 = -1.5;
      a3 = -12.0;
      a4 = 0.1;
      a5 = 1.1;
      a6 = 1.0;
      nu1 = 2.0e-2;
      nu2 = 1e-1;
      
      
      A11 = a1*k - nu1*k^2 + a3;
      A12 = a2*k;
      A21 = a4*k;
      A22 = a5*k - nu2*k^2 + a6;
      
      A = [A11 A12; ...
           A21 A22];

end
%----------------------------------------------------------------------  


