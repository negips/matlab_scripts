%     Plot backflow ratios

clear
clc
close all

load leading_edge_bubble

plot(bubble_time, abs(bubble_ubmax)./bubble_ue, 'LineWidth', 2)
pbaspect([2.5 1 1])
ylabel('$\frac{|U_{b}|}{U_{e}}$', 'FontSize', lafs+4, 'Rotation', 0)
xlabel('$t/T_{osc}$', 'FontSize', lafs)
ylbl = get(gca,'YLabel');
ypos = get(ylbl, 'Position');
set(ylbl,'Position', ypos + [-0.02 0 0])

SaveFig(gcf,'backflow_ratio.eps', 'plots/', 1)
