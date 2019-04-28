% generate the matrix A
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

