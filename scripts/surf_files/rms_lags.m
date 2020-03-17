% Analyzing lag effects in turbulence

clear
clc
close all

mksz = 6;
gr   = (1+sqrt(5))/2;

load rms_all

alpha = ini_aoa + ptch_amp*sin(omega*(time-ptch_start)+ini_phase);


% plot all profiles
nprofs = length(x)

outpos = [0.45 0.05 0.6 0.4];
figure(1)
set(gcf,'Units','Normalized')
set(gcf,'OuterPosition',outpos);
for i=1:nprofs
  X0(:,i) = mean(x{i},2);
  Y0(:,i) = mean(y{i},2);
end

for i=1:nprofs
  plot(X0,Y0, '.b', 'MarkerSize', 4); hold on
end


pr=50;
iy=10;
%plot(x{pr}(iy,1),y{pr}(iy,1), 'or', 'MarkerSize', mksz, 'LineWidth', 2);
plot(X0(iy,pr),Y0(iy,pr), 'or', 'MarkerSize', mksz, 'LineWidth', 2);


var=uu{pr}(iy,:);
outpos2 = [0.65 0.6 0.4 0.4*gr];
figure(2)
set(gcf,'Units','Normalized')
set(gcf,'OuterPosition',outpos2);
plot(alpha,var); hold on

tfine = linspace(time(1),time(end),2000);
vfine = interp1(time,var,tfine);
ar = areafcn(tfine,vfine,0);

lag0 = 0;
[oplag area exitflag output] = fminsearch(@(lag) areafcn(tfine,vfine,lag),lag0);
lag=oplag;

%lag=-1.6;
afine = getalpha(tfine,lag);
plot(afine,vfine,'r.')



%j=0;
%for i=6:nprofs
%   j=j+1;
%   bern(:,j,:) = 0.5*(UU{j}+VV{j}+WW{j}) + P{j};
%   vort(:,j,:) = Wz{j};
%   X(:,j,:)    = x{j};
%   Y(:,j,:)    = y{j};
%   YN(:,j,:)   = yn{j};
%   varuu(:,j,:)= uu{j};
%   Um(:,j,:)   = U{j};
%   Vm(:,j,:)   = V{j};
%
%end
%
% ts = 1;  % time stamp
% 
% b1=bern(:,:,ts);
% vort1=vort(:,:,ts);
% uu1=varuu(:,:,ts);
% um1=Um(:,:,ts);
% vm1=Vm(:,:,ts);
% 
% x1=X(:,:,ts);
% y1=Y(:,:,ts);
% 
% figure(2)
% surf(x1,y1,um1,'EdgeColor','interp')
% %contour(x1,y1,abs(vort1), 50, 'EdgeColor','interp')
% view(2)
% colorbar






