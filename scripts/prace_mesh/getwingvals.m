% random file
% getting some values from the wing

close all
clear
clc

addpath '../'

load('wing_prabal.mat');

Rec = 400*1000;
rho = 1.;
nu = 1/Rec;

yplusmin = 0.64;
yplusmax = 9;
xplus = 18.8;
zplus = 7;
ifxavg = 1;
ifzavg = 1;
lz=0.1;

N = 9;

l1 = length(top);

for i = 1:l1
     xa(i) = top(i).xa;
     ya(i) = top(i).ya;
     ut(i) = top(i).ut;
     tauw(i) = top(i).tauw;
end

l2 = length(bottom);

for i = 1:l1
     xab(i) = bottom(i).xa;
     yab(i) = bottom(i).ya;
     utb(i) = bottom(i).ut;
     tauwb(i) = bottom(i).tauw;
end


%plot(xa,ut)
%hold on

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
display(['===== For x<=0.1 ====='])

%ind = find(xa<=0.1));
%x1 = xa(ind);
%lstar1 = lstar(ind);
%min(lstar1)
%max(lstar1)
%%plot(x1,lstar1)
lstar1 = interp1(xa,lstar,0.2);
if ifzavg
     deltaz = zplus*max(lstar1)*2/avg_gll;
else
     deltaz = zplus*max(lstar1)*2/max_gll;
end
deltaz_gll = deltaz*max_gll/2; 

display(['delta z (Element Size) throughout the domain = ' num2str(deltaz)])
nelz=lz/deltaz;
display(['No. of elements in z (Nelz) = ' num2str(nelz)])
if ifzavg
     display(['Average delta z (gll_pts) throughout the domain = ' num2str(deltaz_gll)])
else
     display(['Max delta z (gll_pts) throughout the domain = ' num2str(deltaz_gll)])
end    

lstar1 = interp1(xa,lstar,0.1);

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

[val indfup] = min(abs(xa - 0.1));
lstf_up = [lstf_up lstar1*ones(1,indfup)]; 
lst_ind = indfup;

xupf = [xupf xa(1:indfup)];
yupf = [yupf ya(1:indfup)];

%% 0.1<x<=0.6
display(['===== For 0.1<x<=0.6 ====='])
[val ind] = min(abs(xa-0.1));
x1 = xa(ind+1:end);
y1 = ya(ind+1:end);
indfup1 = ind+1;

lstar1 = lstar(ind+1:end);
[val ind] = min(abs(x1-0.6));
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
else
     deltax_el_min = xplus*min(lstar2)*2/max_gll;
     deltax_el_max = xplus*max(lstar2)*2/max_gll ;
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

display(['y1 at the x=0.1 = ' num2str(y_st)])
display(['y1 at the x=0.6 = ' num2str(y_end)])

display(['max delta y at x=0.1 = ' num2str(y_st_max)])
display(['max delta y at x=0.6 = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])
display(' ')

els = xplus*lstar2*2/max_gll;

%plot(x2,els)

%% x>0.6
display(['====== For x>0.6 ======'])

[val ind] = min(abs(xa-0.6));
x1 = xa(ind+1:end);
y1 = ya(ind+1:end);

lstar1 = interp1(xab,lstarb,x1);

indfup1 = ind+1;

[val ind2] = min(abs(x1-0.98));                       % Last point is crappy. 0.98 viually chosen.

lstar2 = lstar1(1:ind2);
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

display(['y1 at the x=0.6 = ' num2str(y_st)])
display(['y1 at the x=0.98 = ' num2str(y_end)])

display(['max delta y at x=0.6 = ' num2str(y_st_max)])
display(['max delta y at x=0.98 = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])
display(' ')

els = xplus*lstar2*2/max_gll;

plot(xupf,lstf_up, '--b')

%

%---------------------------------------------------------------------- 

%% Bottom surface
display('Bottom Surface')
display('--------------')

%% x<0.1
display(['===== For x<=0.1 ====='])

lstar1 = interp1(xab,lstarb,0.1);

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

[val indfbot] = min(abs(xab - 0.1));
lstf_bot = [lstf_bot lstar1*ones(1,indfbot)]; 
lst_ind = indfbot;

xbotf = [xbotf xab(1:indfbot)];

%% 0.1<x<=0.6
display(['===== For 0.1<x<=0.6 ====='])
[val ind] = min(abs(xab-0.1));
x1 = xab(ind+1:end);
y1 = yab(ind+1:end);
lstar1 = lstarb(ind+1:end);
indfbot1 = ind+1;

[val ind] = min(abs(x1-0.6));
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
else
     deltax_el_min = xplus*min(lstar2)*2/max_gll;
     deltax_el_max = xplus*max(lstar2)*2/max_gll ;
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

display(['y1 at the x=0.1 = ' num2str(y_st)])
display(['y1 at the x=0.6 = ' num2str(y_end)])

display(['max delta y at x=0.1 = ' num2str(y_st_max)])
display(['max delta y at x=0.6 = ' num2str(y_end_max)])

display(['min delta x = ' num2str(deltax_el_min)])
display(['max delta x = ' num2str(deltax_el_max)])
display(' ')

els = xplus*lstar2*2/max_gll;

%plot(x2,els)

%% x>0.6
display(['====== For x>0.6 ======'])
[val ind] = min(abs(xab-0.6));
x1 = xab(ind+1:end);
y1 = yab(ind+1:end);

lstar1 = lstarb(ind+1:end);

[val ind] = min(abs(x1-0.98));

lstar2 = lstar1(1:ind);
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

display(['y1 at the x=0.6 = ' num2str(y_st)])
display(['y1 at the x=1.0 = ' num2str(y_end)])

display(['max delta y at x=0.6 = ' num2str(y_st_max)])
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

load('data_contour.mat');

maxeta_ratio = 9;
ifxavg = 1;

epsilon2 = abs(0.5*(Dxx + Dyy + Dzz));
eta2 = ((nu^3)./epsilon2).^(1/4);

st_pt = 1.0;
en_pt = max(X(1,:));

[val ind1] = min(abs(X(1,:)-st_pt));

[val ind2] = min(abs(X(1,:)-en_pt));

X2 = X(:,ind1:ind2);
Y2 = Y(:,ind1:ind2);
eta3 = eta2(:,ind1:ind2);
[mineta_x yind] = min(eta3);

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
moveta_min = slidefun(@min,70,mineta_x);
plot(X2(1,:),moveta_min,'m','LineWidth',1);
moveta_s = smooth(moveta_min,140);                % final eta(x)
plot(X2(1,:),moveta_s,'r','LineWidth',1);
legend({'Eta(x)', 'Lower bounds' 'Smoothened Lower bound'}, 'Location', 'Best', 'FontSize', 12);
ylabel('$\eta$', 'Interpreter', 'Latex', 'FontSize', 16)
title('Kolmogorov scale', 'Interpreter', 'Latex', 'FontSize', 16)

linfit = fit(transpose(X2(1,500:3000)), moveta_s(500:3000), 'poly1');
plot(linfit, 'k');

%% Isotropy

dxx2 = Dxx(:,ind1:ind2);
dyy2 = Dyy(:,ind1:ind2);

xlocs = [];
xlocs = linspace(1.5,5,8);
xlocs = [1.1 xlocs];

figure
hold on
cmap = lines(length(xlocs));
for i = 1:length(xlocs)
     [val ind3] = min(abs(X2(1,:)-xlocs(i)));
     h_dr(i) = plot(Y2(:,ind3),dxx2(:,ind3),'Color', cmap(i,:));
     plot(Y2(:,ind3),dyy2(:,ind3), '--', 'Color',cmap(i,:));
     legs{i} = ['x=' num2str(xlocs(i))];
end
xlim([-0.25 0.5]);
legend(h_dr,legs);
title('Dxx, Dyy');

%% Taking deltax ~ 2*deltay
%% (dx*dy*dz)^(1/3) = 10*eta; ==> dx = [(10*eta)^(1/3)/dz*2]^(1/2);
%% ==> el_s = dx/2*avg_gll;
h_x = maxeta_ratio*moveta_s;

% Using volume criterion
%deltax = [(h_x.^3)/deltaz_gll*2].^(1/2);           % delta z calculated from the wing surface

% using length criterion
deltax = h_x;

if ifxavg
     el_s_x = deltax*2/avg_gll;
else
     el_s_x = deltax*2/max_gll;
end

display('Wake')
display('-----------')

val1 = interp1(X2(1,:),el_s_x,2.0);
display(['delta x @x=2.0 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,3.0);
display(['delta x @x=3.0 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,4.0);
display(['delta x @x=4.0 == ' num2str(val1)])

val1 = interp1(X2(1,:),el_s_x,5.0);
display(['delta x @x=5.0 == ' num2str(val1)])

figure
plot(X2(1,:),el_s_x);
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('Element Sizes', 'Interpreter', 'Latex', 'FontSize', 16);


%% Kolmogorov scale above the airfoil
st_pt = 0.5;
en_pt = 1.2;

[val ind1] = min(abs(X(1,:)-st_pt));

[val ind2] = min(abs(X(1,:)-en_pt));

X3 = X(:,ind1:ind2);
Y3 = Y(:,ind1:ind2);
eta4 = eta2(:,ind1:ind2);

st_pt = 0;
en_pt = 0.15;

[val ind1] = min(abs(Y3(:,1)-st_pt));
[val ind2] = min(abs(Y3(:,1)-en_pt));
X4 = X3(ind1:ind2,:);
Y4 = Y3(ind1:ind2,:);
eta5 = eta4(ind1:ind2,:);
h_xa = maxeta_ratio*eta5;

% Some queries @

[val yind] = min(abs(Y4(:,1)-0.05));

xpts = [0.6 0.7 0.8 0.9 1.0];
cmap = lines(length(xpts));

figure
hold on
legs =[];
for i=1:length(xpts);
     [val ind1] = min(abs(X4(1,:)-xpts(i)));
     yvals = Y4(1:yind,ind1);
     hvals = h_xa(1:yind,ind1);
     plot(yvals,hvals,'Color', cmap(i,:));
     legs{i} = ['x = ' num2str(xpts(i))];

end

title('"h" above the airfoil')
legend(legs)


