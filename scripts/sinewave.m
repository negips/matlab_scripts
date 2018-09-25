%% Script to make figure for phase angle.

clear
clc
%close all

fsi = importdata('fsi_io.out');
l1 = length(fsi.data(:,1));
iskip = 100;
ind2 = find(fsi.data(:,2)>80,1);
ind = 1:iskip:l1;

figure(1)
set(gcf,'Color', 'none')
set(gcf,'renderer','painters', 'Units','normalized');
set(gcf,'OuterPosition', [0.25 0.25 0.45 0.40])
%ax1=axes;
plot(fsi.data(ind,2),fsi.data(ind,4), 'm', 'LineWidth', 4);
set(gca,'XAxisLocation','origin')
set(gca,'XTick', []);
set(gca,'YTick', []);
hold on
%set(gca, 'XTick', [0 pi/2 pi 3*pi/2 2*pi])
%set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '2\pi'}, 'TickLabelInterpreter', 'latex', 'FontSize', 60)
%set(gca,'YTick', [])
%set(gca,'Color','none')
%hold on
%xlim([0,2*pi+0.1])
%ylim([-5.35 5.35])

%set(gca,'Color','none')
%set(gcf,'Color', 'none')
filename = ['sinewave2.eps']; 
SaveFig(gcf,filename,'plots/',1)

