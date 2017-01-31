%    Compare NEK GLL point ditribution
%    Get dx+ values


clear
clc
close all

addpath '../'

%load('wing_prabal.mat');

Rec = 100*1000;
nu = 1/Rec;

yplusmin = 0.64;
yplusmax = 12;
xplus = 18.8;
zplus = 9;
ifxavg = 1;
ifzavg = 1;

N = 11;

xfoil = importdata('integral_vals_re100k.out');
l1 = length(xfoil.data(:,7));
x0 = xfoil.data(1,2);
cnt=0;
rho=1;
while(x0~=0)
     cnt=cnt+1;
     xa(cnt) = xfoil.data(cnt,2);
     ya(cnt) = xfoil.data(cnt,3);
     tauw(cnt) = abs(xfoil.data(cnt,7))*0.5;
     ut(cnt) = sqrt(tauw(cnt)/rho);

     x0 = xfoil.data(cnt+1,2);
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
lstf_up = ut/nu;

xab = transpose(xfoil.data(cnt+1:end,2));
yab = transpose(xfoil.data(cnt+1:end,3));
tauwb = transpose(abs(xfoil.data(cnt+1:end,7))*0.5);
utb = sqrt(tauwb/rho);
lstf_bot = utb/nu;

%%
nek_dat = importdata('pitch.N11');

% Divide into top and bottom halves.
% Using a simplistic criterion. Maybe something more robut later...

[nx ny] = size(nek_dat.data);

up_cnt =0;
lo_cnt =0;
co_upp = [];
co_low = [];
for ii = 1:nx

     if (nek_dat.data(ii,5)>=0)
%          if up_cnt==0 || (abs(nek_dat.data(ii,1)-co_upp(up_cnt,1))>1e-7 || abs(nek_dat.data(ii,2)-co_upp(up_cnt,2))>1e-7)
               up_cnt=up_cnt+1;
               co_upp(up_cnt,:) = [nek_dat.data(ii,1) nek_dat.data(ii,2)];
%          end
     else
%          if lo_cnt==0 || (abs(nek_dat.data(ii,1)-co_low(lo_cnt,1))>1e-7 || abs(nek_dat.data(ii,2)-co_low(lo_cnt,2))>1e-7)
               lo_cnt=lo_cnt+1;
               co_low(lo_cnt,:) = [nek_dat.data(ii,1) nek_dat.data(ii,2)];
%          end
     end
end


[co_upp(:,1) I] = sort(co_upp(:,1));
co_upp(:,2) = co_upp(I,2);

[co_low(:,1) I] = sort(co_low(:,1));
co_low(:,2) = co_low(I,2);

plot(co_upp(:,1),co_upp(:,2))
hold on
plot(co_low(:,1),co_low(:,2), 'r')


%% Distance along the curves
dup = diff(co_upp);
dlo = diff(co_low);

dpts = [];
nx = length(xa);
for i = 1:nx
     dpts(i) = sqrt((xa(i)-0)^2 + (ya(i)-0)^2);
end

dpts_neku = [];
[nx ny] = size(co_upp);
dr_nek = [];
for i = 1:nx
     dpts_neku(i) = sqrt((co_upp(i,1)-0)^2 + (co_upp(i,2)-0)^2);
     if i~=nx
          dr_neku(i) = sqrt((co_upp(i+1,1)-co_upp(i,1))^2 + (co_upp(i+1,2)-co_upp(i,2))^2);
     end
end

dpts_nekl = [];
[nx ny] = size(co_low);
for i = 1:nx
     dpts_nekl(i) = sqrt((co_low(i,1)-0)^2 + (co_low(i,2)-0)^2);
     if i~=nx
          dr_nekl(i) = sqrt((co_low(i+1,1)-co_low(i,1))^2 + (co_low(i+1,2)-co_low(i,2))^2);
     end
end


lstar_gllup = interp1(dpts,lstf_up,dpts_neku);
lstar_gll = lstar_gllup(1:end-1);
eta_r_u = dr_neku./lstar_gll;
figure
plot(co_upp(21:end-1,1),eta_r_u(21:end), '.');      % Some weird values before this
hold on

lstar_gllbot = interp1(dpts,lstf_bot,dpts_nekl);
lstar_gll = lstar_gllbot(1:end-1);
eta_r_l = dr_nekl./lstar_gll;
plot(co_low(1:end-1,1),eta_r_l(1:end), 'r.');
ylim([0 35])
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('\Delta x^{+}', 'Interpreter', 'tex')
legend({'Top surface', 'Bottom surface'}, 'Location', 'Best');

