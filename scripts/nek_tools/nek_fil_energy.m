% Change in energy due to nek filtering

clear
clc
close all

lafs=20;


alphas=1.0;
N=10;
Un = 1;
Un_2 = 1;

En = 1/(2*N+1);
En_2 = 1/(2*(N-2)+1);

%Un_range = linspace(-100,100,25);
%Un_2_range = linspace(-100,100,25);
%
%[Un Un_2] = meshgrid(Un_range,Un_2_range);
%
%%dE = (Un_2.^2 - (Un_2 - alpha*Un).^2).*En_2 + (Un.^2 - ((1-alpha).*Un).^2).*En;

%surf(Un,Un_2,dE,'EdgeColor', 'none', 'FaceColor', 'interp')
%xlabel('U_{N}')
%ylabel('U_{N-2}')

ratio = linspace(-5,5,500);         % ratio (Un-2/Un)
alphas = linspace(1e-6,1.00,500);

[a r] = meshgrid(alphas,ratio);

dE_sign = ((r.^2 - (r+a).^2).*En_2 + (1 - (1-a).^2).*En);
dE_sign2 = sign(dE_sign);

s = surf(a,r,dE_sign2,'EdgeColor', 'none', 'FaceColor', 'interp');
cmap = flipud(colormap('gray'));
colormap(cmap)
caxis([-1 3])
%colorbar
hold on
%alpha(s,0.5)
%contour(a,r,dE_sign, [0 0], 'LineColor', 'k', 'LineWidth', 2)
view(2)
%line([0 max(alpha)], [0 0], [10 10], 'LineStyle', '--', 'LineWidth', 1.5)
%ylabel('$\hat{U}_{N-2}/\hat{U}_{N}$', 'Interpreter', 'Latex')
ylabel('$r$', 'FontSize', lafs)
xlabel('$\alpha$', 'FontSize', lafs)

dim = [0.4 0.5 0.20 0.3];
str = 'Energy Source';
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
set(an1,'FontSize',16)
% set(an1, 'EdgeColor', [1 1 1])
set(an1,'Margin', 10)

dim = [0.4 0.10 0.20 0.3];
str = 'Energy Sink';
an2 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
set(an2,'FontSize',16)
set(an2,'Margin', 10)
% set(an2, 'LineWidth', 0)


SaveFig(gcf, 'filter_delta_energy.eps', 'plots/', 1)


