%% Plot integral quantities

clear
clc
close all

lafs=16;

area=0.25*1;
U0=1.0;
rho=1.0;
norm = 0.5*rho*U0^2*area;

k=0.5;
Chord=1.0;
semichord=Chord/2;
omega=k*U0/semichord;
Tosc=2*pi/omega;

cl=importdata('cl.out');

ind = find(cl.data(:,1)==0);
cl.data(ind,:) = [];

ind = (cl.data(:,2)<=17.75);  % 17.87
t1 = cl.data(find(ind),2);
cl1 = cl.data(find(ind),4);
ind = (cl.data(:,2)>18.85);
t2 = cl.data(find(ind),2);
cl2 = cl.data(find(ind),4);

hack = load('hack.mat');
fac1=7.6621;
fac2=7.9972;
ind1 = hack.time_hist>17.87;
ind2 = hack.time_hist<18.85;
tmid = hack.time_hist(find(ind1.*ind2));
clmid = hack.cl_hist(find(ind1.*ind2));
fac=linspace(fac1,fac2,length(clmid));
clscaled = clmid.*(fac')*norm;

timeall=[t1; tmid; t2];
clall  =[cl1; clscaled; cl2];

ax1=axes;
plot(timeall/Tosc,clall/norm,'b', 'LineWidth', 2,'Parent', ax1);
set(ax1, 'YLim', [1.1 1.6])
set(ax1, 'Xlim', [min(cl.data(:,2)) max(cl.data(:,2))]/Tosc)
ylabel('$C_{L}$', 'FontSize', lafs, 'Parent', ax1)
xlabel('$t/T_{osc}$', 'FontSize', lafs, 'Parent', ax1)
%set(ax1,'YColor', [0 0 1])
lg1 = legend({'$C_{L}$'}, 'FontSize', lafs, 'Location', 'SouthWest', 'Color', 'none', 'Box', 'on');

ax2=axes;
phi = omega*timeall;
alpha = 6.7 + 1.3*sin(phi);
plot(timeall/Tosc,alpha, '--k', 'LineWidth', 1.5, 'Parent', ax2)
%set(ax2,'Xlim',[1 25]);
set(ax2,'PlotBoxAspectRatio', [2.5 1 1])
set(ax1,'PlotBoxAspectRatio', [2.5 1 1])
set(ax2,'XAxisLocation', 'top');
set(ax2,'YAxisLocation', 'right');
set(ax2,'Color', 'none')
set(ax2,'XColor', [1 1 1])
%set(ax2,'YColor', [1 0 0])
set(ax2, 'Xlim', [min(cl.data(:,2)) max(cl.data(:,2))]/Tosc)
set(ax2,'XTick', [])
ylabel('$\alpha[^{o}]$', 'FontSize', lafs, 'Parent', ax2)
lg2 = legend({'$\alpha$'}, 'FontSize', lafs, 'Location', 'SouthEast', 'Color', 'none', 'Box', 'on');
set(ax2,'YLim', [4 9])

SaveFig(gcf,'cl-time-alpha.eps', 'plots/',1)

figure(2)
ax3=axes;
ind=find(timeall>4.0*pi);
cl2=clall(ind);
alpha2=alpha(ind);
plot(alpha2,cl2/norm, 'LineWidth', 2, 'Parent',ax3)
ph = arrowh(alpha2,cl2/norm,'k',[200,100],[10 30 40 60 80]);
xlabel('$\alpha[^{o}]$', 'FontSize', lafs, 'Parent', ax3)
ylabel('$C_{L}$', 'FontSize', lafs, 'Parent', ax3)
SaveFig(gcf,'cl-alpha.eps', 'plots/',1)



% figure(3)
% ax3=axes;
% ind=find(cl.data(:,2)>4.0*pi);
% cl2=cl.data(ind,4);
% OMEGA = 1.3*omega*cos(phi);
% OMEGA2=OMEGA(ind);
% plot(OMEGA2,cl2/norm, 'LineWidth', 1.5, 'Parent',ax3)
% xlabel('\Omega', 'FontSize', 20, 'Parent', ax3)
% ylabel('C_{L}', 'FontSize', 20, 'Parent', ax3)
% %SaveFig(gcf,'cl-omega.eps', 'plots/',1)

