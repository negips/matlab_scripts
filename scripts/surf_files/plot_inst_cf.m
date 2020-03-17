% Plot the instantaneous Cf

clear
clc
close all

load('re750k_surfaceN9.mat');
axfs=30;
lafs=44;

%destn = 'plots_test/';
destn = '~/workstation/git_kth/forced_pitching/paper/imgs2/';

tnorm = (surf_t9(:,1)-ptch_start)/Tosc;

gr = (1+sqrt(5))/2;
figpos = [0.25 0.25 0.35 0.35*gr];

% Plot close to min tr point
figure(1);
set(gcf,'Units','normalized')
set(gcf,'OuterPosition',figpos)
j=2408;
plot(surf_x9(j,:),surf_v9(j,:), 'LineWidth', 3)
ylim([-1.5 5]*10^-3)
set(gca,'FontSize',axfs)
set(gca,'LineWidth',1.5)
xlabel('$x/c$', 'FontSize',lafs)
ylabel('$\tau_{w}$', 'FontSize',lafs)
SaveFig(gcf,'pitch750k_t05_70',destn,1)
t0 = tnorm(j)


% Plot close to min tr point
figure(2);
set(gcf,'Units','normalized')
set(gcf,'OuterPosition',figpos)
j=4410;
plot(surf_x9(j,:),surf_v9(j,:), 'LineWidth', 3)
ylim([-1.5 5]*10^-3)
set(gca,'FontSize',axfs)
set(gca,'LineWidth',1.5)

xlabel('$x/c$', 'FontSize',lafs)
ylabel('$\tau_{w}$', 'FontSize',lafs)
SaveFig(gcf,'pitch750k_t06_21',destn,1)
t0 = tnorm(j)


figure(3)
plot(tnorm)
