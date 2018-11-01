%% Creating curved side data (mid point)

clear
clc
close all

c1 = 'k';
c2 = 'k';
c3 = [0.4 0.4 0.4]; % gray
ls1 = '-';
ls2 = '--';
ls3 = ':';
lw1 = 2;
lw2 = 2;
lw3 = 0.5;

lafs = 20;

th = 0.02; % half thickness
rad = 0.5;

theta = asin(th/rad);

s = linspace(theta,2*pi-theta,1000);

xc=rad*cos(s);
yc=rad*sin(s);
figure(1);
set(gcf,'Units','Normalized')
set(gcf,'OuterPosition', [0.25 0.25 0.55 0.45])

plot(xc,yc, ls1, 'Color', c1, 'LineWidth', lw1); hold on
%xlim([-0.6 0.6])
%ylim([-0.6 0.6])

xst = rad*cos(theta);
yst = th;

% splitter plate
xs1 = [xst xst+1];
ys1 = [yst yst];

xs2 = [xst+1 xst+1];
ys2 = [th -th];

xs3 = [xst+1 xst];
ys3 = [-yst -yst];


xsplitter = [xs1 xs2 xs3];
ysplitter = [ys1 ys2 ys3];

plot(xsplitter,ysplitter, ls1, 'Color', c1, 'LineWidth', lw1)

% Axis
xax = [-1 2];
yax = [0 0];
plot(xax,yax, ls3, 'Color', c3, 'LineWidth', lw3)


% Rotate
angle = 15*pi/180;

Rot = [cos(angle)  sin(angle);
       -sin(angle) cos(angle)];

xy2 = Rot*[xc; yc];
xc2 = xy2(1,:);
yc2 = xy2(2,:);

xy2 = Rot*[xsplitter; ysplitter];

xsplitter2 = xy2(1,:);
ysplitter2 = xy2(2,:);

plot(xc2,yc2, ls2, 'Color', c2, 'LineWidth', lw2)
plot(xsplitter2,ysplitter2, ls2, 'Color', c2, 'LineWidth', lw2)

% Axis
xy2 = Rot*[xax; yax];
xax2 = xy2(1,:);
yax2 = xy2(2,:);
plot(xax2,yax2, ls3, 'Color', c3, 'LineWidth', lw3)

xlabel('$x$', 'FontSize', lafs)
ylabel('$y$', 'FontSize', lafs)
ylim([-0.62 0.62])
%axis equal

% Annotation
andim = [0.27 0.52 0.1 0.1];
ann = annotation('textbox',andim,'String','$\theta$');
set(ann, 'LineStyle', 'none')
set(ann, 'Interpreter', 'latex')
set(ann, 'FontSize', 16)

SaveFig(gcf,'cyl_rot.eps', 'plots/', 1)


