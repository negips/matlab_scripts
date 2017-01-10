% Stability of BDF2/EXT2 scheme

clear
clc
close all

x1 = [-2:0.05:2];
y1 = [-2:0.05:2]*1i;

for i = 1:length(x1)

x = x1(i);
for j = 1:length(y1)
y = y1(j);
phi = x + y;

a = 3.00001;
b = 4*(phi-1);
c = 1-2*phi;

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

%h1 = figure;
%plot(real_r,imag_r, 'o')
%ylabel('Imag(G)')
%xlabel('Real(G)')

[X,Y] = meshgrid(x1,y1/1i);
 
h2 = figure;
surf(X,Y,abs(transpose(r1)),'EdgeColor', 'none')
hold on
contour(X,Y,abs(transpose(r1)),[1,1],'LineWidth',4,'LineColor','k','ShowText', 'on')
xlabel('Real(\phi)')
ylabel('Imag(\phi)')
colorbar
view([0 -90])
%xlim([-0.1 1.5])
%ylim([-1 1])
% max(max(abs(r1)))

hold on
h3 = figure;
surf(X,Y,abs(transpose(r2)),'EdgeColor', 'none')
hold on
contour(X,Y,abs(transpose(r2)),[1,1],'LineWidth',4,'LineColor','k','ShowText', 'on')
xlabel('Real(\phi)')
ylabel('Imag(\phi)')
%grid on
colorbar
view([0 -90])
%max(max(abs(r2)))



