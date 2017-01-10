% Stability of BDF2/EXT2 scheme

clear
clc
close all

beta_arr = [2:0.1:5]*64*10^-4;
l = length(beta_arr);

gamma_arr = [-1:0.05:0];
k1 = length(gamma_arr);

pe = 100000;
k = 1;
C = 1;

for i = 1:l
beta = beta_arr(i);

theta = 90/180*3.1415;

for j = 1:k1

gamma = gamma_arr(j);

a = 3.00001;
b = (2i*beta*sin(theta) -gamma -4);
c = (1 - 1i*beta*sin(theta) +0*gamma);

rts = roots([a b c]);
r1(i,j) = rts(1);
r2(i,j) = rts(2);

end

end

real_r = real([r1 r2]);

[l1 l2] = size(real_r);

real_r = reshape(real_r,l1*l2,1);

imag_r = imag([r1 r2]);
imag_r = reshape(imag_r,l1*l2,1);

h1 = figure;
plot(real_r,imag_r, 'o')
ylabel('Imag(G)')
xlabel('Real(G)')

[X,Y] = meshgrid(beta_arr,gamma_arr);
 
h2 = figure;
surf(X,Y,abs(transpose(r1)))
xlabel('beta')
ylabel('gamma')
colorbar
max(max(abs(r1)))

hold on
h3 = figure;
surf(X,Y,abs(transpose(r2)))
xlabel('beta')
ylabel('gamma')
grid on
colorbar
max(max(abs(r2)))



