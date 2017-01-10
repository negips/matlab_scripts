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

a = 2.0000;
b = 3*(6*phi+1);
c = -6*(3*phi+1);
d = 1+6*phi;

rts = roots([a b c d]);
r1(i,j) = rts(1);
r2(i,j) = rts(2);
r3(i,j) = rts(3);

end

end

real_r = real([r1 r2 r3]);

[l1 l2] = size(real_r);

real_r = reshape(real_r,l1*l2,1);

imag_r = imag([r1 r2 r3]);
imag_r = reshape(imag_r,l1*l2,1);

%h1 = figure;
%plot(real_r,imag_r, 'o')
%ylabel('Imag(G)')
%xlabel('Real(G)')

[X,Y] = meshgrid(x1,y1/1i);
 
h2 = figure;
surf(X,Y,abs(transpose(r1)),'EdgeColor', 'none')
hold on
%contour(X,Y,abs(transpose(r1)),[1 1],'LineWidth',2,'LineColor','k','ShowText', 'on')
hold on
xlabel('Real(\phi)')
ylabel('Imag(\phi)')
colorbar
%view([0 -90])
title('R1');
%xlim([-0.1 1.5])
%ylim([-1 1])
% max(max(abs(r1)))
print(h2,'-depsc', '-r300','h2.eps')

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
title('R2')

h4 = figure;
surf(X,Y,abs(transpose(r3)),'EdgeColor', 'none')
hold on
contour(X,Y,abs(transpose(r3)),[1,1],'LineWidth',4,'LineColor','k','ShowText', 'on')
xlabel('Real(\phi)')
ylabel('Imag(\phi)')
title('R3')
%grid on
colorbar
view([0 90])


