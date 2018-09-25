% random file
% getting some values from the wing

close all
clear
clc

addpath '../'

alpha = 0/180*pi;
uy = 1.*sin(alpha);
ux = 1*cos(alpha);
uinf = 1.;

%surfdata = importdata('ed36F128+14_aoa66.bl-upper');
%l1=length(surfdata.data(:,2));
%xa = fliplr(transpose(surfdata.data(:,2)));
%ya = fliplr(transpose(surfdata.data(:,3)));
%cf = fliplr(transpose(surfdata.data(:,7)));
%ut = sqrt(abs(cf)/2*uinf^2);
%
%surfdata = importdata('ed36F128+14_aoa66.bl-lower');
%l1=length(surfdata.data(:,2));
%xab = transpose(surfdata.data(:,2));
%yab = transpose(surfdata.data(:,3));
%cfb = transpose(surfdata.data(:,7));
%utb = sqrt(abs(cfb)/2*uinf^2);
%
%ind = find(xab>1);
%xab(ind) = [];
%yab(ind) = [];
%cfb(ind) = [];
%utb(ind) = [];

Rec = 77*1000;
rho = 1.0;
nu = 1/Rec;
nu_nek=nu;

yplusmin = 0.64;
yplusmax = 12;
xplus = 18.0;
zplus = 12;
ifxavg = 1;
ifzavg = 1;
lz = 0.2;

% Data from xfoil

%xfoil = importdata('integral_vals_re100k.out');
xfoil = importdata('integral_vals_re77k_naca_alpha50.out');
l1 = length(xfoil.data(:,7));
x0 = xfoil.data(1,2);
cnt=0;
dx0 = diff(xfoil.data(:,2));
while(dx0(cnt+1)<0)
  cnt=cnt+1;
  xa(cnt) = xfoil.data(cnt,2);
  ya(cnt) = xfoil.data(cnt,3);
  tauw(cnt) = abs(xfoil.data(cnt,7))*0.5;
  ut(cnt) = sqrt(tauw(cnt)/rho);

%  x0 = xfoil.data(cnt+1,2);
end

cnt=cnt+1;
xa(cnt) = xfoil.data(cnt,2);
ya(cnt) = xfoil.data(cnt,3);
tauw(cnt) = abs(xfoil.data(cnt,7))*0.5;
ut(cnt) = sqrt(tauw(cnt)/rho);

xa = fliplr(xa);
ya = fliplr(ya);
tauw = fliplr(tauw);
ut = fliplr(ut);

xab = transpose(xfoil.data(cnt+1:end,2));
yab = transpose(xfoil.data(cnt+1:end,3));
tauwb = transpose(abs(xfoil.data(cnt+1:end,7))*0.5);
utb = sqrt(tauwb/rho);


N = 11;

figure
plot(xa,ut)
hold on
plot(xab,utb,'-r')
ylabel('$u_{\tau}$', 'Interpreter', 'Latex')
xlabel('$x$', 'Interpreter', 'Latex')
legend({'Top', 'Bottom'}, 'Location', 'Best')

lstf_up = [];
lstf_bot = [];
xupf = [];
xbotf = [];
yupf = [];
ybotf = [];

lstar = nu./ut;
lstarb = nu./utb;

figure
plot(xa,lstar)
hold on
plot(xab,lstarb,'r')
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('$l^{*}$', 'Interpreter', 'Latex', 'FontSize', 16);
legend({'Top', 'Bottom'}, 'Interpreter', 'Latex', 'FontSize', 16, 'Location', 'Best');

gll = lglnodes(N);
gll_diff = abs(diff(gll));

min_gll = min(gll_diff);
max_gll = max(gll_diff);
avg_gll = mean(gll_diff);

%% Top surface
display('Top Surface')
display('-----------')

%% x<0.1
xref=0.185;

display(['===== For x<=' num2str(xref) ' ====='])

%ind = find(xa<=0.1));
%x1 = xa(ind);
%lstar1 = lstar(ind);
%min(lstar1)
%max(lstar1)
%%plot(x1,lstar1)
lstar1 = interp1(xa,lstar,xref);             % Has to be selected manually
if ifzavg
     deltaz = zplus*max(lstar1)*2/avg_gll;
else
     deltaz = zplus*max(lstar1)*2/max_gll;
end
nelz = lz/deltaz;
deltaz_gll = deltaz*max_gll/2; 

display(['delta z (Element Size) throughout the domain = ' num2str(deltaz)])
display(['No of Elements in Z = ' num2str(nelz) ' for Lz =' num2str(lz)])

if ifzavg
     display(['Average delta z (gll_pts) throughout the domain = ' num2str(deltaz_gll)])
else
     display(['Max delta z (gll_pts) throughout the domain = ' num2str(deltaz_gll)])
end    

lstar1 = interp1(xa,lstar,xref);

min_el_size = yplusmin*min(lstar1)*2/min_gll;
max_el_size = yplusmax*max(lstar1)*2/max_gll;

if ifxavg==1
     deltax = xplus*max(lstar1)*2/avg_gll;
else
     deltax = xplus*max(lstar1)*2/max_gll;
end

y1 = min_el_size;
ymax = max_el_size;

display(['y1 at the wall = ' num2str(y1)])
display(['max delta y = ' num2str(ymax)])
display(['delta x = ' num2str(deltax)])
display (' ')

[val indfup] = min(abs(xa - xref));
lstf_up = [lstf_up lstar1*ones(1,indfup)]; 
lst_ind = indfup;

xupf = [xupf xa(1:indfup)];
yupf = [yupf ya(1:indfup)];

%% 0.1<x<=0.6
xref1=xref;
xref2=0.60;
display(['===== For ' num2str(xref1) '<x<=' num2str(xref2) ' ====='])
[val ind] = min(abs(xa-xref1));
x1 = xa(ind+1:end);
y1 = ya(ind+1:end);
indfup1 = ind+1;

lstar1 = lstar(ind+1:end);
[val ind] = min(abs(x1-xref2));
x2 = x1(1:ind);
y2 = y1(1:ind);
lstar2 = lstar1(1:ind);
indfup2 = lst_ind + ind;

lstf_up = [lstf_up lstar2];
lst_ind=indfup2;
xupf = [xupf x2];
yupf = [yupf y2];

if ifxavg==1
     deltax_el_min = xplus*min(lstar2)*2/avg_gll;
     deltax_el_max = xplus*max(lstar2)*2/avg_gll ;

     deltax_el_1 = xplus*lstar2(1)*2/avg_gll;
     deltax_el_2 = xplus*lstar2(end)*2/avg_gll ;

else
     deltax_el_min = xplus*min(lstar2)*2/max_gll;
     deltax_el_max = xplus*max(lstar2)*2/max_gll ;

     deltax_el_1 = xplus*lstar2(1)*2/max_gll;
     deltax_el_2 = xplus*lstar2(end)*2/max_gll ;

end

y1 = yplusmin*min(lstar2)*2/min_gll;
y2 = yplusmin*max(lstar2)*2/min_gll;

y_st = yplusmin*lstar2(1)*2/min_gll;
y_st_max = yplusmax*lstar2(1)*2/min_gll;

y_end = yplusmin*lstar2(end)*2/min_gll;
y_end_max = yplusmax*lstar2(end)*2/min_gll;

ymax = yplusmax*max(lstar2)*2/max_gll;

display(['y1 at the wall(min) = ' num2str(y1)])
%display(['y1 at the wall(max) = ' num2str(y2)])

display(['y1 at the x=' num2str(xref1) ' = ' num2str(y_st)])
display(['y1 at the x=' num2str(xref2) ' = ' num2str(y_end)])

display(['max delta y at x=' num2str(xref1) ' = ' num2str(y_st_max)])
display(['max delta y at x=' num2str(xref2) ' = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])

display(['delta x @ x= ' num2str(xref1) ' ==> ' num2str(deltax_el_min)])
display(['delta x @ x= ' num2str(xref2) ' ==> ' num2str(deltax_el_max)])

display(' ')

els = xplus*lstar2*2/max_gll;

%plot(x2,els)

%% x>0.6
xref3=xref2;
display(['====== For x>' num2str(xref3) ' ======'])

[val ind] = min(abs(xa-xref3));
x1 = xa(ind+1:end);
y1 = ya(ind+1:end);

lstar1 = lstar(ind+1:end);

indfup1 = ind+1;

ind2=length(x1);                                      % Constant grid spacing after 0.7 

lstar2 = lstar1(1:ind2);
lstar2 = min(lstar2)*ones(1,ind2);
x2 = x1(1:ind2);
y2 = y1(1:ind2);

indfup2 = lst_ind + ind2;
lstf_up = [lstf_up lstar2];
lst_ind=indfup2;
xupf = [xupf x2];
yupf = [yupf y2];

y1 = yplusmin*min(lstar2)*2/min_gll;
y2 = yplusmin*max(lstar2)*2/min_gll;

y_st = yplusmin*lstar2(1)*2/min_gll;
y_st_max = yplusmax*lstar2(1)*2/min_gll;

y_end = yplusmin*lstar2(end)*2/min_gll;
y_end_max = yplusmax*lstar2(end)*2/min_gll;

ymax = yplusmax*max(lstar2)*2/max_gll;

if ifxavg==1
     deltax_el_min = xplus*min(lstar2)*2/avg_gll;
     deltax_el_max = xplus*max(lstar2)*2/avg_gll ;
else
     deltax_el_min = xplus*min(lstar2)*2/max_gll;
     deltax_el_max = xplus*max(lstar2)*2/max_gll ;
end

display(['y1 at the wall(min) = ' num2str(y1)])
%display(['y1 at the wall(max) = ' num2str(y2)])

display(['y1 at the x=' num2str(xref) ' = ' num2str(y_st)])
display(['y1 at the x=0.98 = ' num2str(y_end)])

display(['max delta y at x=' num2str(xref) ' = ' num2str(y_st_max)])
display(['max delta y at x=0.98 = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])
display(' ')

els = xplus*lstar2*2/max_gll;

plot(xupf,lstf_up, '--b')

%---------------------------------------------------------------------- 

%% Bottom surface
display('Bottom Surface')
display('--------------')

%% x<0.1
xref=0.1325;
display(['===== For x<=' num2str(xref) ' ====='])

lstar1 = interp1(xab,lstarb,xref);

min_el_size = yplusmin*min(lstar1)*2/min_gll;
max_el_size = yplusmax*max(lstar1)*2/max_gll;

if ifxavg==1
     deltax = xplus*max(lstar1)*2/avg_gll;
else
     deltax = xplus*max(lstar1)*2/max_gll;
end

y1 = min_el_size;
ymax = max_el_size;

display(['y1 at the wall = ' num2str(y1)])
display(['max delta y = ' num2str(ymax)])
display(['delta x = ' num2str(deltax)])
display (' ')

[val indfbot] = min(abs(xab - xref));
lstf_bot = [lstf_bot lstar1*ones(1,indfbot)]; 
lst_ind = indfbot;

xbotf = [xbotf xab(1:indfbot)];

%% 0.1<x<=0.6
xref1=0.1325;
xref2=0.7;
display(['===== For ' num2str(xref1) '<x<=' num2str(xref2) ' ====='])
[val ind] = min(abs(xab-xref1));
x1 = xab(ind+1:end);
y1 = yab(ind+1:end);
lstar1 = lstarb(ind+1:end);
indfbot1 = ind+1;

[val ind] = min(abs(x1-xref2));
x2 = x1(1:ind);
y2 = y1(1:end);
lstar2 = lstar1(1:ind);

indfbot2 = lst_ind + ind;

lstf_bot = [lstf_bot lstar2];
lst_ind=indfbot2;
xbotf = [xbotf x2];
ybotf = [ybotf y2];

plot(xbotf,lstf_bot, '--r');

if ifxavg==1
     deltax_el_min = xplus*min(lstar2)*2/avg_gll;
     deltax_el_max = xplus*max(lstar2)*2/avg_gll ;

     deltax_el_1 = xplus*lstar2(1)*2/avg_gll;
     deltax_el_2 = xplus*lstar2(end)*2/avg_gll ;

else
     deltax_el_min = xplus*min(lstar2)*2/max_gll;
     deltax_el_max = xplus*max(lstar2)*2/max_gll ;

     deltax_el_1 = xplus*lstar2(1)*2/max_gll;
     deltax_el_2 = xplus*lstar2(end)*2/max_gll ;

end

y1 = yplusmin*min(lstar2)*2/min_gll;
y2 = yplusmin*max(lstar2)*2/min_gll;

y_st = yplusmin*lstar2(1)*2/min_gll;
y_st_max = yplusmax*lstar2(1)*2/min_gll;

y_end = yplusmin*lstar2(end)*2/min_gll;
y_end_max = yplusmax*lstar2(end)*2/min_gll;

ymax = yplusmax*max(lstar2)*2/max_gll;

display(['y1 at the wall(min) = ' num2str(y1)])
%display(['y1 at the wall(max) = ' num2str(y2)])

display(['y1 at the x=' num2str(xref1) ' = ' num2str(y_st)])
display(['y1 at the x=' num2str(xref2) ' = ' num2str(y_end)])

display(['max delta y at x=' num2str(xref1) ' = ' num2str(y_st_max)])
display(['max delta y at x=' num2str(xref2) ' = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])

display(['delta x @ x= ' num2str(xref1) ' ==> ' num2str(deltax_el_min)])
display(['delta x @ x= ' num2str(xref2) ' ==> ' num2str(deltax_el_max)])

display(' ')

els = xplus*lstar2*2/max_gll;

%plot(x2,els)

%% x>0.6
xref=0.7;
display(['====== For x>' num2str(xref) ' ======'])
[val ind] = min(abs(xab-xref));
x1 = xab(ind+1:end);
y1 = yab(ind+1:end);

lstar1 = lstarb(ind+1:end);

ind = length(x1);

lstar2 = lstar1(1:ind);
lstar2 = min(lstar2)*ones(1,ind);
x2 = x1(1:ind);
y2 = y1(1:ind);

indfbot2 = lst_ind + ind2;
lstf_bot = [lstf_bot lstar2];
lst_ind=indfbot2;
xbotf = [xbotf x2];
ybotf = [ybotf y2];

plot(xbotf,lstf_bot, '--r')

y1 = yplusmin*min(lstar2)*2/min_gll;
y2 = yplusmin*max(lstar2)*2/min_gll;

y_st = yplusmin*lstar2(1)*2/min_gll;
y_st_max = yplusmax*lstar2(1)*2/min_gll;

y_end = yplusmin*lstar2(end)*2/min_gll;
y_end_max = yplusmax*lstar2(end)*2/min_gll;

ymax = yplusmax*max(lstar2)*2/max_gll;

if ifxavg==1
     deltax_el_min = xplus*min(lstar2)*2/avg_gll;
     deltax_el_max = xplus*max(lstar2)*2/avg_gll ;
else
     deltax_el_min = xplus*min(lstar2)*2/max_gll;
     deltax_el_max = xplus*max(lstar2)*2/max_gll ;
end

display(['y1 at the wall(min) = ' num2str(y1)])
%display(['y1 at the wall(max) = ' num2str(y2)])

display(['y1 at the x=' num2str(xref) ' = ' num2str(y_st)])
display(['y1 at the x=1.0 = ' num2str(y_end)])

display(['max delta y at x=' num2str(xref) ' = ' num2str(y_st_max)])
display(['max delta y at x=1.0 = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])
display(' ')

els = xplus*lstar2*2/max_gll;

save('lstar_wing.mat', 'xbotf', 'ybotf', 'xupf', 'yupf', 'lstf_bot', 'lstf_up')

%figure
%plot(x2,els)

%======================================================================


% Wake regions values
display('Wake')
display('-----------')

% load('../prace_mesh/data_contour.mat');

%rans = importdata('epsilon-rsm-03.dat');
rans = importdata('naca_77k_eps.dat');

filsgs_eta = 9;
ifxavg = 1;
ifsgs = 1;
sgs_ramp_start=1.25;
sgs_full=2.00;
sgs_eta=20;


%epsilon2 = abs(0.5*(Dxx + Dyy + Dzz));
X = rans.data(:,2);
Y = rans.data(:,3);
epsilon = rans.data(:,4);
nu_rans = (1.7894e-5)/1.225;
uinf_rans = 1.1248;
epsilon_norm = epsilon*(nu_rans/(uinf_rans^4));
epsilon_nek = epsilon_norm/(nu/1^4);
%eta = (nu_nek^(3/4))/(epsilon_nek.^(1/4));
eta=((nu_nek^3)./epsilon_nek).^(1/4);

st_pt = 1.0;
en_pt = 5.2;

intp_pts = 1000;
xvec = linspace(st_pt,en_pt,intp_pts);

st_pt = -2;
en_pt = 2;
yvec = linspace(st_pt,en_pt,intp_pts);

[X2 Y2] = meshgrid(xvec,yvec);
eta2 = griddata(X,Y,eta,X2,Y2);

[mineta_x yind] = min(eta2);

ylocs_m=zeros(1,length(yind));
for i=1:length(ylocs_m)
  ylocs_m(i) = Y2(yind(i),i);
end
%scatter(X2(1,:),ylocs_m,mineta_x);
figure
plot(X2(1,:),ylocs_m,'-*', 'MarkerSize', 12)
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel('$y$', 'Interpreter', 'Latex', 'FontSize', 16)
title('Location of peak dissipation in the wake', 'Interpreter', 'Latex', 'FontSize', 16)


figure
plot(X2(1,:),mineta_x);
hold on
% Get bounds
% Since from the rans it is already smooth enough. I don't need smoothening, unlike the data from DNS.
% Hence smoothing/moving minima interval set to 1.
interval = 1;
moveta_min = slidefun(@min,interval,mineta_x);
plot(X2(1,:),moveta_min,'m','LineWidth',1);

interval=1;
moveta_s = smooth(moveta_min,interval);                % final eta(x) used.
plot(X2(1,:),moveta_s,'r','LineWidth',1);
ylabel('$\eta$')

legend({'Eta(x)', 'Lower bounds' 'Smoothened Lower bound'}, 'Location', 'Best', 'FontSize', 12);
ylabel('$\eta$', 'Interpreter', 'Latex', 'FontSize', 16)
title('Kolmogorov scale', 'Interpreter', 'Latex', 'FontSize', 16)

linfit = fit(transpose(X2(1,:)), moveta_s, 'poly1');
plot(linfit, 'k');

%% Isotropy

%-------------------- 

if (ifsgs)
  smoothvals = smoothstep(X2(1,:),sgs_ramp_start,sgs_full);
  maxeta_ratio = filsgs_eta + smoothvals.*(sgs_eta-filsgs_eta);
  figure(20)
  plot(X2(1,:),maxeta_ratio)
else
  maxeta_ratio=filsgs_eta;
end
h_x = transpose(maxeta_ratio).*moveta_s;


% Using volume criterion
% deltax = [(h_x.^3)/deltaz_gll*2].^(1/2);           % delta z calculated from the wing surface

% using length criterion
deltax = h_x;

if ifxavg
  el_s_x = deltax*2/avg_gll;
else
  el_s_x = deltax*2/max_gll;
end

val1 = interp1(X2(1,:),el_s_x,1.25);
display(['delta x @x=1.25 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,2.0);
display(['delta x @x=2.0 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,3.0);
display(['delta x @x=3.0 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,4.0);
display(['delta x @x=4.0 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,5.2);
display(['delta x @x=5.2 == ' num2str(val1)])

figure
plot(X2(1,:),el_s_x);
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('Element Sizes', 'Interpreter', 'Latex', 'FontSize', 16);

dh_el = el_s_x;
dh_ratio = dh_el(2:end)./dh_el(1:end-1);
figure(21)
plot(X2(1,2:end),dh_ratio)
grid on
grid minor

%% Kolmogorov scale above the airfoil

% st_pt = -0.2;
% en_pt = 1.2;
% 
% intp_pts = 200;
% xvec = linspace(st_pt,en_pt,intp_pts);
% 
% st_pt = 0;
% en_pt = 0.5;
% yvec = linspace(st_pt,en_pt,intp_pts);
% 
% [X3 Y3] = meshgrid(xvec,yvec);
% eta3 = griddata(X,Y,eta,X3,Y3);
% 
% % Bit retarded this way, but it was adapted from the code that used DNS data.
% % You know, whatever works.
% st_pt = 0;
% en_pt = 0.4;
% 
% [val ind1] = min(abs(Y3(:,1)-st_pt));
% [val ind2] = min(abs(Y3(:,1)-en_pt));
% X4 = X3(ind1:ind2,:);
% Y4 = Y3(ind1:ind2,:);
% eta4 = eta3(ind1:ind2,:);
% h_xa = maxeta_ratio.*eta4;
% 
% % Some queries @
% 
% [val yind] = min(abs(Y4(:,1)-0.05));
% 
% xpts = [0.6 0.7 0.8 0.9 1.0];
% cmap = lines(length(xpts));
% 
% figure
% hold on
% legs =[];
% for i=1:length(xpts);
%      [val ind1] = min(abs(X4(1,:)-xpts(i)));
%      yvals = Y4(1:yind,ind1);
%      hvals = h_xa(1:yind,ind1);
%      plot(yvals,hvals,'Color', cmap(i,:));
%      legs{i} = ['x = ' num2str(xpts(i))];
% end
% 
% title('"h" above the airfoil')
% legend(legs)


