% Reynolds stresses

clear
clc
close all

filename='xy.stat';
[stat,x,y,Re,fltype,dstar,rlam,spanv,wAstat,sumw]=readxystats(filename);

[X,Y]=meshgrid(x,y);

x_start = 1000;
x_end   = 11000;
y_start = 0;
y_end   = 250;
ind1=x>=x_start;
ind2=x<=x_end;
indx = find(ind1.*ind2);
x2 = x(indx);
ind1=y>=y_start;
ind2=y<=y_end;
indy = find(ind1.*ind2);
y2 = y(indy);

[X2,Y2]=meshgrid(x2,y2);

total_p = 0.5*(stat.u2 + stat.v2 + stat.w2) + stat.p;

bl_criteria=0.99;
delta99 = [];
for i=1:length(x)
  delta99(i)=interp1(stat.u(:,i),y,bl_criteria);
  delta995(i)=interp1(stat.u(:,i),y,0.995);
  Ue(i)=interp1(y,stat.u(:,i),delta99(i));
  Ve(i)=interp1(y,stat.v(:,i),delta99(i));
end

dtmp = delta99(indx);
for i=1:length(indx)
  Y3(:,i) = Y2(:,i)/dtmp(i);
end

figure(100)
surf(X2,Y3,stat.urms(indy,indx),'EdgeColor', 'none', 'FaceColor', 'interp');
view(2)
ylim([0 1.2]);
colorbar
hold on

%figure(100)
%plot3(x(indx),delta99(indx),100*ones(size(delta99(indx))), 'LineWidth', 2, 'Color', 'k' )
%view(2)


ymin = min(y);
ymax = max(y);
Ly = ymax - ymin;

cheb_u_deficit = chebfun(1.0 - stat.u);
dstar = sum(cheb_u_deficit)*Ly;
cheb_theta_deficit = chebfun((1.0 - stat.u).*stat.u);
theta = sum(cheb_theta_deficit)*Ly;

figure(1)
plot(x,delta99)
hold on
plot(x,delta995, '--b')
plot(x,dstar, 'r')
plot(x,theta, 'k')
legend({'\delta_{99}', '\delta_{995}', '\delta_{*}', '\theta' }, 'FontSize', 14, 'Location', 'NorthWest')

figure(2)
plot(x,dstar./theta)
ylabel('H')
xlim([1000 12000])


figure(2)
plot(x,dstar./delta99)
hold on
plot(x,theta./delta99, 'r')
%plot(x,delta99./delta995, 'k')
xlim([1000 11000])


% generate wall normal profiles

xin = [1000 2000 5000 8000]
ifynorm = 1;









