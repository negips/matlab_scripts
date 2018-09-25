%% Script to make figure for phase angle.

clear
clc
close all

t=linspace(0,2*pi+0.5,1000);
omega = 1.0;
phi = omega*t;
alpha_0 = 6.7;
dalpha = 1.3;

alpha = alpha_0 + dalpha*sin(phi);

figure(1)
set(gcf,'Color', 'none')
set(gcf,'renderer','painters', 'Units','normalized');
set(gcf,'OuterPosition', [0.25 0.25 0.45 0.40])
ax1=axes;
plot(phi,alpha-alpha_0, 'LineWidth', 4);
set(gca,'XAxisLocation','origin')
set(gca, 'XTick', [0 pi/2 pi 3*pi/2 2*pi])
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '2\pi'}, 'TickLabelInterpreter', 'latex', 'FontSize', 60)
set(gca,'YTick', [])
set(gca,'Color','none')
hold on
xlim([0,2*pi+0.1])
ylim([-1.35 1.35])


positions = [0.09 0.2 0.3 0.47 0.55 0.63 0.82 1.01];
npos = length(positions);

dim=[0.6 0.45 0.4 0.4];

for i=1:npos
  ph = 2*pi*positions(i);
  pos = dalpha*sin(ph);
  pl = plot(ph,pos, 'or','LineWidth',2, 'MarkerSize',30, 'MarkerFaceColor', 'r');
  str=['$\alpha=', num2str(pos+alpha_0,3) '^{o}$']; 
  str=['$', num2str(pos+alpha_0,3) '^{o}$']; 
  an = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'Interpreter', 'latex', 'FontSize', 100, 'EdgeColor', 'none');
  set(gca,'Color','none')
  set(gcf,'Color', 'none')
  filename = ['phase_angle',num2str(100*positions(i),3), '.eps']; 
  SaveFig(gcf,filename,'plots/',1)
  delete(pl)
  delete(an)
end 


