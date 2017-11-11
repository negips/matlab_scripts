%% Plot integral quantities

clear
clc
close all

area=0.15*1;
U0=1.0;
rho=1.0;
norm = 0.5*rho*U0^2*area;

ax1=axes;
CLimp=importdata('cl.out');

ind = find(CLimp.data(:,1)==0);
CLimp.data(ind,:) = [];
ind1 = CLimp.data(:,2)>10.0;
ind2 = CLimp.data(:,2)<30.0;
ind=find(ind1.*ind2);

CLimp.data = CLimp.data(ind,:);
time = CLimp.data(:,2);

cl = CLimp.data(:,4)/norm;

plot(time,cl,'b', 'LineWidth', 1.5,'Parent', ax1);
%set(ax1, 'YLim', [1.1 1.6])
set(ax1, 'Xlim', [min(time) max(time)])
ylabel('C_{L}', 'FontSize', 20, 'Parent', ax1)
xlabel('Time', 'FontSize', 20, 'Parent', ax1)
%set(ax1,'YColor', [0 0 1])

%ax2=axes;
k=0.4;
U0=1.0;
Chord=1.0;
semichord=Chord/2;
omega=k*U0/semichord;
phase_shift=-pi/2;
ptch_start=6.0;
phi = omega*(time-ptch_start);
alpha = 3.4 + 1.0*sin(phi+phase_shift);

%plot(CLimp.data(:,2),alpha, '--k', 'LineWidth', 1.5, 'Parent', ax2)
%%set(ax2,'Xlim',[1 25]);
%set(ax2,'XAxisLocation', 'top');
%set(ax2,'YAxisLocation', 'right');
%set(ax2,'Color', 'none')
%set(ax2,'XColor', [1 1 1])
%%set(ax2,'YColor', [1 0 0])
%set(ax2, 'Xlim', [min(CLimp.data(:,2)) max(CLimp.data(:,2))])
%set(ax2,'XTick', [])
%ylabel('\alpha^{o}', 'FontSize', 20, 'Parent', ax2)

%SaveFig(gcf,'CLimp-time-alpha.eps', 'plots/',1)

figure(2)
ax3=axes;
%ind=find(CLimp.data(:,2)>4.0*pi);
ind=[1:length(alpha)];
time2=time(ind);
cl2=cl(ind);
alpha2=alpha(ind);
plot(alpha2,cl2, 'LineWidth', 1.5, 'Parent',ax3)
xlabel('\alpha^{o}', 'FontSize', 20, 'Parent', ax3)
ylabel('C_{L}', 'FontSize', 20, 'Parent', ax3)
SaveFig(gcf,'750k_cl-alpha_sim.eps', 'plots/',1)
hold on


%figure(3)
%ax3=axes;
%ind=find(cl.data(:,2)>4.0*pi);
%cl2=cl.data(ind,4);
%OMEGA = 1.3*omega*cos(phi);
%OMEGA2=OMEGA(ind);
%plot(OMEGA2,cl2/norm, 'LineWidth', 1.5, 'Parent',ax3)
%xlabel('\Omega', 'FontSize', 20, 'Parent', ax3)
%ylabel('C_{L}', 'FontSize', 20, 'Parent', ax3)
%SaveFig(gcf,'cl-omega.eps', 'plots/',1)


%% Fitting model

xfoil = importdata('test_ed36f128+14_re7.5e5.dat');
xfoilalpha = xfoil.data(:,1);
xfoilcl = xfoil.data(:,2);

%expd = load('14_static_models_765k.mat');
%
%xfoilalpha = expd.alpha-0.2;
%
%cz_24 = interp1(expd.alpha,expd.cz,2.4);
%xfoilcl = expd.cz + cl2(1)-cz_24;


figure(2)
plot(xfoilalpha,xfoilcl, '--k')

theta=0;
amp2=1;
par0(1)=-20;            % phase lag (quasi-steady)
par0(2)=0.1;            % integration
par0(3)=-50;            % phase-lag added mass
par0(4)=0.05;           % cl offset

cz=cl2;
alpha0 = 3.4;
dalpha = 1.0;
omega = omega; 

options=optimset('MaxFunEvals',1000,'MaxIter',10000,'TolX',1e-8,'Tolfun',1e-8);
[par,fval,exitflag,output] = fminsearch(@(par) PhaseLag(par,time2,cz,omega,alpha0,dalpha,xfoilalpha,xfoilcl), par0,options);
%[par,fval,exitflag,output] = fminsearch(@(par) PhaseLag(par,time2,cz,omega,alpha0,dalpha,xfoilalpha,xfoilcl), par0);


phi=par(1)*pi/180;
intg_const=par(2);
theta = par(3)*pi/180;
ofst = par(4); 

phase_lag = -pi/2;

pitch = dalpha*sin(omega*(time2-6.0) + phi + phase_lag);
alpha_lagg = alpha0 + pitch;
added_mass = abs(intg_const)*cos(omega*(time2-6.0) + theta + phase_lag);

cl_lagg = interp1(xfoilalpha,xfoilcl,alpha_lagg,'linear');
cl_pred = ofst + added_mass + cl_lagg;

figure(2)
plot(alpha,cl_pred, 'r')





