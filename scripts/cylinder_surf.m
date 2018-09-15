% make a visualization for a cylinder surface


clear
clc
close all

theta = linspace(0,2*pi,5000);

r = 1.;

x = r*cos(theta);
y = r*sin(theta);


figure(1)
plot(x,y, 'k', 'LineWidth', 4); hold on

% x2 = [x, fliplr(x)];
% inBetween = [y, fliplr(curve2)];
% fill(x2, inBetween, 'g');
%fill(x,y, [0.75 0.75 0.75])
f1 = fill(x,y, [0. 0. 0.]);

x2 = x;
y2 = y+0.01;
plot(x2,y2, '--k', 'LineWidth', 2)

x3 = [x2 fliplr(x)];
y3 = [y2 fliplr(y)];
f2 = fill(x3,y3, [0.85 0.85 0.85], 'LineStyle', 'none');


xlim([-0.15 0.15])
ylim([0.97 1.05])

xbox = [0.25 0.35];
ybox = [0.7 0.45];
ann1 = annotation('textarrow',xbox,ybox,'String','$\Gamma$', 'Interpreter', 'latex');
set(ann1,'FontSize', 20)

xbox = [0.60 0.55];
ybox = [0.75 0.56];
ann2 = annotation('textarrow',xbox,ybox,'String','Perturbed surface position');
set(ann2,'FontSize', 16)

%set(gca,'Visible', 'off')
ylabel('$y$', 'FontSize', 24)
xlabel('$x$', 'FontSize', 24)

SaveFig(gcf,'cylinder_surf', 'plots/', 1)


