%     Test Boostconv and spsnap

clear
clc
close all

%rng('default');

n = 50;    % Matrix Size

% No of eigenvalues greater than 1.
np = 2;
% How unstable do we want the eigenvalues?
uns_range = 0.1;

theta = pi/2*rand(np,1);
e1  = (cos(theta) + uns_range*rand(np,1)) + 1i*sin(theta);
uns_e = [e1; conj(e1)];
% Set first eigenvector = 1 signifying beseflow

ne2 = n-length(e1);
st_range = 0.5;
theta = 2*pi*rand(ne2,1);
damp = 0.1 + 0.85*rand(ne2,1);
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
plot(cos(theta),sin(theta), '--r')
xlabel('$\lambda_{r}$', 'FontSize', 24)
ylabel('$\lambda_{i}$', 'FontSize', 24)


bkryl = 15;

x0 = rand(neig,1) + 1i*rand(neig);

