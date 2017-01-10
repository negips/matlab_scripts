%    Compare NEK GLL point ditribution
%    Get dx+ values


clear
clc
close all

addpath '../'

%load('wing_prabal.mat');

Rec = 400*1000;
nu = 1/Rec;

yplusmin = 0.64;
yplusmax = 12;
xplus = 18.8;
zplus = 7;
ifxavg = 1;
ifzavg = 1;

N = 11;

%l1 = length(top);
%
%for i = 1:l1
%     xa(i) = top(i).xa;
%     ya(i) = top(i).ya;
%     ut(i) = top(i).ut;
%     tauw(i) = top(i).tauw;
%end
%
%l2 = length(bottom);
%
%for i = 1:l1
%     xab(i) = bottom(i).xa;
%     yab(i) = bottom(i).ya;
%     utb(i) = bottom(i).ut;
%     tauwb(i) = bottom(i).tauw;
%end

load('lstar_wing.mat');

xa = xupf;
ya = yupf;
xab = xbotf;
yab = ybotf;

nek_dat = importdata('surf_data.dat');

% Divide into top and bottom halves.
% Using a simplistic criterion. Maybe something more robust later...

[nx ny] = size(nek_dat.data);

up_cnt =0;
lo_cnt =0;
co_upp = [];
co_low = [];
for ii = 1:nx

     if (nek_dat.data(ii,2)>=0)
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

