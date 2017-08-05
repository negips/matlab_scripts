% Create streamlines from the bl data

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

figure(100)
%surf(X2,Y2,stat.u(indy,indx),'EdgeColor', 'none', 'FaceColor', 'interp');
%surf(X2,Y2,atan(stat.v(indy,indx)./stat.u(indy,indx))*180/pi,'EdgeColor', 'none', 'FaceColor', 'interp');
surf(X2,Y2,total_p(indy,indx),'EdgeColor', 'none', 'FaceColor', 'interp');
view(2)
colorbar
hold on

bl_criteria=0.99;
delta99 = [];
for i=1:length(x)
  delta99(i)=interp1(stat.u(:,i),y,bl_criteria);
  delta995(i)=interp1(stat.u(:,i),y,0.995);
  Ue(i)=interp1(y,stat.u(:,i),delta99(i));
  Ve(i)=interp1(y,stat.v(:,i),delta99(i));
  
  ttp = total_p(:,i);
  ttp_max = max(ttp);
  ttp = ttp + (1-ttp_max);
  delta99_p(i)=interp1(ttp,y,0.99);
end

figure(100)
plot3(x(indx),delta99(indx),100*ones(size(delta99(indx))), 'LineWidth', 2, 'Color', 'k' )
view(2)


ymin = min(y);
ymax = max(y);
Ly = ymax - ymin;

cheb_u_deficit = chebfun(1.0 - stat.u);
dstar = sum(cheb_u_deficit)*Ly;
cheb_theta_deficit = chebfun((1.0 - stat.u).*stat.u);
theta = sum(cheb_theta_deficit)*Ly;

%cheby = 2*y/Ly -1;
%for i=1:length(cheby)
%  clc
%  i
%  dstar_prof(i,:)=sum(cheb_u_deficit,-1,cheby(i));
%end

figure(1)
plot(x,delta99)
hold on
plot(x,delta995, '--b')
plot(x,delta99_p, ':b')
plot(x,dstar, 'r')
plot(x,theta, 'k')
legend({'\delta_{99}', '\delta_{995}', '\delta_{99}^{p}','\delta_{*}', '\theta' }, 'FontSize', 14, 'Location', 'NorthWest')

figure(2)
plot(x,dstar./theta)
xlim([1000 12000])

%x_seed = 1000; 
%nseed = 5;
%dr = 1e-2;
%maxverts=828000;
%blh_xin = interp1(x,delta99,x_seed);
%y_seed = linspace(1e-4,1.0*blh_xin,nseed);
%
%x_sl = [];
%y_sl = [];
%y_sl2 = [];
%xsl_max = 12000;
%ysl_max = 400;
%
%for iseed = 1:nseed
%  xsl=x_seed;
%  ysl=y_seed(iseed);
%  xsl_arr = [xsl];
%  ysl_arr = [ysl];
%  npts=0;    
%  while ((xsl<xsl_max) && (ysl<ysl_max))
%    clc  
%    npts=npts+1;
%    iseed  
%    [npts xsl]
%    [xsl ysl] = sl_rk4(X,Y,stat.u,stat.v,dr,xsl,ysl);
%    xsl_arr = [xsl_arr; xsl];
%    ysl_arr = [ysl_arr; ysl];
%    if (mod(npts,1000)==0)
%      figure(20)
%      plot(xsl_arr,ysl_arr)
%      pause(0.1)
%    end    
%  end
%  slines(iseed).x=xsl_arr;
%  slines(iseed).y=ysl_arr;
%  slines(iseed).xw = x;
%  slines(iseed).delta99 = delta99;        
%end     
%
%save('slines.mat', 'slines')

% figure(100)
% for i=1:nseed
%   tmp=find(ind1(:,i).*ind2(:,i));
%   plot3(x_sl(tmp,i),y_sl(tmp,i),100*ones(size(x_sl(tmp,i))), 'LineWidth', 1.5, 'Color', 'k','LineStyle', '--' )
% end
% view(2)

% figure(4)
% plot(x_sl,y_sl2)
% title('Normalized Streamlines')

%% Manually calculated streamlines








