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

Rec = 77*1000;
rho = 1.0;
nu = 1/Rec;
nu_nek=nu;

yplusmin = 0.3;
yplusmax = 8;
xplus = 10.0;
zplus = 8;
ifxavg = 1;
ifzavg = 1;
lz = 0.10;

% Data from xfoil

xfoil = importdata('integral_vals_re77k_naca_alpha00.out');
l1 = length(xfoil.data(:,7));
x0 = xfoil.data(1,2);
cnt=0;

[val ind0] = min(abs(xfoil.data(:,2))-0.);
indt = 1:ind0;
indb = ind0+1:length(xfoil.data(:,2));

xa = xfoil.data(indt,2)';
ya = xfoil.data(indt,3)';
tauw = transpose(abs(xfoil.data(indt,7))*0.5);
ut = sqrt(tauw/rho);
xa = fliplr(xa);
ya = fliplr(ya);
tauw = fliplr(tauw);
ut = fliplr(ut);

xab = xfoil.data(indb,2)';
yab = xfoil.data(indb,3)';
tauwb = transpose(abs(xfoil.data(indb,7))*0.5);
utb = sqrt(tauwb/rho);


N = 10;

figure(1)
plot(xa,ut)
hold on
plot(xab,utb,'-r')
ylabel('$u_{\tau}$', 'Interpreter', 'Latex')
xlabel('$x$', 'Interpreter', 'Latex')
legend({'Top', 'Bottom'}, 'Location', 'Best')


top_breaks = [0.18 0.56];
topsections = length(top_breaks) + 1;
top_breaks = [0 top_breaks 1.01];      % Just add end points
top_min    = [0 0 1];        % If take the maximum value with this section
top_xref   = [0.1778 0 0];   % If use this reference point for calculation of grid

lstf_up = [];
xupf = [];
yupf = [];

lstar = nu./ut;

figure(2)
lst = plot(xa,lstar); hold on

gll = lglnodes(N);
gll_diff = abs(diff(gll));

min_gll = min(gll_diff);
max_gll = max(gll_diff);
avg_gll = mean(gll_diff);

%% Top surface
display('Top Surface')
display('-----------')

for ii = 1:topsections

  disp(['Section: ' num2str(top_breaks(ii)) '<= x < ', num2str(top_breaks(ii+1))])
  if (top_min(ii))
    disp('Taking the maximum value of Tauw')
  end  
  if (top_xref(ii)>0)
    disp(['Taking Tauw from point: x=', num2str(top_xref(ii))])
  end  
  disp('----------------------------------------')

  ind1 = xa>=top_breaks(ii);
  ind2 = xa<top_breaks(ii+1);
  ind3 = find(ind1.*ind2);
  x1 = xa(ind3);
  y1 = ya(ind3);

  if (top_min(ii))                        % Take the minimum value
    lstar1 = interp1(xa,lstar,x1);
    lstar_min = min(lstar1);
    lstar1 = zeros(size(x1)) + lstar_min;
  elseif (top_xref(ii)>0)                 % Take from reference point
    lstar_ref = interp1(xa,lstar,top_xref(ii));
    lstar1 = zeros(size(x1)) + lstar_ref;
  else
    lstar1 = interp1(xa,lstar,x1);
  end

  min_el_size = yplusmin*min(lstar1)*2/min_gll;
  max_el_size = yplusmax*max(lstar1)*2/max_gll;
  
  if ifxavg==1
    deltax = xplus*max(lstar1)*2/avg_gll;
  else
    deltax = xplus*max(lstar1)*2/max_gll;
  end

  if ifxavg==1
    deltax_el_min = xplus*min(lstar1)*2/avg_gll;
    deltax_el_max = xplus*max(lstar1)*2/avg_gll ;
  
    deltax_el_1 = xplus*lstar1(1)*2/avg_gll;
    deltax_el_2 = xplus*lstar1(end)*2/avg_gll ;
  else
    deltax_el_min = xplus*min(lstar1)*2/max_gll;
    deltax_el_max = xplus*max(lstar1)*2/max_gll ;
  
    deltax_el_1 = xplus*lstar1(1)*2/max_gll;
    deltax_el_2 = xplus*lstar1(end)*2/max_gll ;
  end

%  ymin = min_el_size;
%  ymax = max_el_size;
  
  yp1 = yplusmin*min(lstar1)*2/min_gll;
  yp2 = yplusmin*max(lstar1)*2/min_gll;
  
  y_st = yplusmin*lstar1(1)*2/min_gll;
  y_st_max = yplusmax*lstar1(1)*2/min_gll;
  
  y_end = yplusmin*lstar1(end)*2/min_gll;
  y_end_max = yplusmax*lstar1(end)*2/min_gll;
  
  ymax = yplusmax*max(lstar1)*2/max_gll;
  
  display(['y^{+} at the wall(min) = ' num2str(yp1)])
  display(['y^{+} at the wall(max) = ' num2str(yp2)])
  
  display(['y^{+} at the x=' num2str(top_breaks(ii)) ' = ' num2str(y_st)])
  display(['y^{+} at the x=' num2str(top_breaks(ii+1)) ' = ' num2str(y_end)])
  
%  display(['max delta y at x=' num2str(top_breaks(ii)) ' = ' num2str(y_st_max)])
%  display(['max delta y at x=' num2str(top_breaks(ii+1)) ' = ' num2str(y_end_max)])
  
  display(['min delta x = ' num2str(deltax_el_min)])
  display(['max delta x = ' num2str(deltax_el_max)])
  
  display(['delta x @ x= ' num2str(top_breaks(ii)) ' ==> ' num2str(deltax_el_min)])
  display(['delta x @ x= ' num2str(top_breaks(ii+1)) ' ==> ' num2str(deltax_el_max)])
  display (' ')

  lstf_up = [lstf_up lstar1]; 
  
  xupf = [xupf x1];
  yupf = [yupf y1];

end   % ii=1:topsections

pla = plot(xupf,lstf_up, '--b');
%---------------------------------------------------------------------- 
%---------------------------------------------------------------------- 

bot_breaks = [0.18 0.56];
botsections = length(bot_breaks) + 1;
bot_breaks = [0 bot_breaks 1];      % Just add end points
bot_min    = [0 0 1];        % If take the maximum value with this section
bot_xref   = [0.1778 0 0];   % If use this reference point for calculation of grid

lstf_bot = [];
xbotf = [];
ybotf = [];

lstarb = nu./utb;

plot(xab,lstarb,'r')

%% Bottom surface
display('Bottom Surface')
display('-----------')

for ii = 1:botsections

  disp(['Section: ' num2str(bot_breaks(ii)) '<= x < ', num2str(bot_breaks(ii+1))])
  if (bot_min(ii))
    disp('Taking the maximum value of Tauw')
  end  
  if (bot_xref(ii)>0)
    disp(['Taking Tauw from point: x=', num2str(bot_xref(ii))])
  end  
  disp('----------------------------------------')

  ind1 = xab>=bot_breaks(ii);
  ind2 = xab<bot_breaks(ii+1);
  ind3 = find(ind1.*ind2);
  x1 = xab(ind3);
  y1 = yab(ind3);

  if (bot_min(ii))                        % Take the minimum value
    lstar1 = interp1(xab,lstarb,x1);
    lstar_min = min(lstar1);
    lstar1 = zeros(size(x1)) + lstar_min;
  elseif (bot_xref(ii)>0)                 % Take from reference point
    lstar_ref = interp1(xab,lstarb,bot_xref(ii));
    lstar1 = zeros(size(x1)) + lstar_ref;
  else
    lstar1 = interp1(xab,lstarb,x1);
  end
  
  min_el_size = yplusmin*min(lstar1)*2/min_gll;
  max_el_size = yplusmax*max(lstar1)*2/max_gll;
  
  if ifxavg==1
    deltax = xplus*max(lstar1)*2/avg_gll;
  else
    deltax = xplus*max(lstar1)*2/max_gll;
  end

  if ifxavg==1
    deltax_el_min = xplus*min(lstar1)*2/avg_gll;
    deltax_el_max = xplus*max(lstar1)*2/avg_gll ;
  
    deltax_el_1 = xplus*lstar1(1)*2/avg_gll;
    deltax_el_2 = xplus*lstar1(end)*2/avg_gll ;
  else
    deltax_el_min = xplus*min(lstar1)*2/max_gll;
    deltax_el_max = xplus*max(lstar1)*2/max_gll ;
  
    deltax_el_1 = xplus*lstar1(1)*2/max_gll;
    deltax_el_2 = xplus*lstar1(end)*2/max_gll ;
  end
  
  yp1 = yplusmin*min(lstar1)*2/min_gll;
  yp2 = yplusmin*max(lstar1)*2/min_gll;
  
  y_st = yplusmin*lstar1(1)*2/min_gll;
  y_st_max = yplusmax*lstar1(1)*2/min_gll;
  
  y_end = yplusmin*lstar1(end)*2/min_gll;
  y_end_max = yplusmax*lstar1(end)*2/min_gll;
  
  ymax = yplusmax*max(lstar1)*2/max_gll;
  
  display(['y^{+} at the wall(min) = ' num2str(yp1)])
  display(['y^{+} at the wall(max) = ' num2str(yp2)])
  
  display(['y^{+} at the x=' num2str(bot_breaks(ii)) ' = ' num2str(y_st)])
  display(['y^{+} at the x=' num2str(bot_breaks(ii+1)) ' = ' num2str(y_end)])
  
  display(['max delta y at x=' num2str(bot_breaks(ii)) ' = ' num2str(y_st_max)])
  display(['max delta y at x=' num2str(bot_breaks(ii+1)) ' = ' num2str(y_end_max)])
  
  display(['min delta x = ' num2str(deltax_el_min)])
  display(['max delta x = ' num2str(deltax_el_max)])
  
  display(['delta x @ x= ' num2str(bot_breaks(ii)) ' ==> ' num2str(deltax_el_min)])
  display(['delta x @ x= ' num2str(bot_breaks(ii+1)) ' ==> ' num2str(deltax_el_max)])
  display (' ')

  lstf_bot = [lstf_bot lstar1]; 
  
  xbotf = [xbotf x1];
  ybotf = [ybotf y1];

end   % ii=1:botsections

plb = plot(xbotf,lstf_bot, '--r');

xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('$l^{*}$', 'Interpreter', 'Latex', 'FontSize', 16);
legend([pla plb], {'Top', 'Bottom'}, 'Interpreter', 'Latex', 'FontSize', 16, 'Location', 'Best');

% Find smallest resolution for Delta z
lst_all = [lstf_up lstf_bot];
if ifzavg
  deltaz = zplus*min(lst_all)*2/avg_gll;
else
  deltaz = zplus*min(lst_all)*2/max_gll;
end
nelz = lz/deltaz;
deltaz_gll = deltaz*min_gll/2; 

disp('----------------------------------------')
display(['delta z (Element Size) throughout the domain = ' num2str(deltaz)])
display(['No of Elements in Z is ' num2str(nelz) ' for Lz =' num2str(lz)])

if ifzavg
  display(['Average delta z (gll_pts) throughout the domain = ' num2str(deltaz_gll)])
else
  display(['Max delta z (gll_pts) throughout the domain = ' num2str(deltaz_gll)])
end    

save('lstar_wing.mat', 'xbotf', 'ybotf', 'xupf', 'yupf', 'lstf_bot', 'lstf_up')

%return
%figure
%plot(x2,els)

%======================================================================

% break

% Wake regions values
display('Wake')
display('-----------')

% load('../prace_mesh/data_contour.mat');

%rans = importdata('epsilon-rsm-03.dat');
rans = importdata('naca_150k_eps.dat');

filsgs_eta = 9;
ifxavg = 1;
ifsgs = 1;
sgs_ramp_start=1.25;
sgs_full=2.00;
sgs_eta=9;


%epsilon2 = abs(0.5*(Dxx + Dyy + Dzz));
X = rans.data(:,2);
Y = rans.data(:,3);
epsilon = rans.data(:,4);
nu_rans = (1.7894e-5)/1.225;
uinf_rans = 2.27925;
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


