% Analyzing lag effects in turbulence

clear
clc
close all

mksz = 6;
gr   = (1+sqrt(5))/2;

ts = 54;  % time stamp
pr=63;
iy=60;
ifsave = 0;

load rms_all

alpha = ini_aoa + ptch_amp*sin(omega*(time-ptch_start)+ini_phase);

% plot all profiles
nprofs = length(x)

outpos = [0.35 0.05 0.4*gr 0.6];
figure(1)
set(gcf,'Units','Normalized')
set(gcf,'OuterPosition',outpos);
for i=1:nprofs
  X0(:,i) = mean(x{i},2);
  Y0(:,i) = mean(y{i},2);
end

Z0 = 10*ones(size(X0));

for i=1:nprofs
%  plot3(X0(1:75,:),Y0(1:75,:), Z0(1:75,:), '.b', 'MarkerSize', 3); hold on
  plot3(x{i}(1:75,ts),y{i}(1:75,ts), Z0(1:75,1),'.k', 'MarkerSize', 3); hold on
end
plot3(x{pr}(iy,ts),y{pr}(iy,ts), Z0(iy,pr), 'om', 'MarkerSize', mksz, 'LineWidth', 2);
%plot3(X0(iy,pr),Y0(iy,pr), Z0(iy,pr), 'or', 'MarkerSize', mksz, 'LineWidth', 2);
ylabel('$y$')
xlabel('$x$')
ttl = ['Profile No:',num2str(pr), '; iy:', num2str(iy)];
title(ttl);

foil = importdata('ed36F128+14.dat');
fx = foil.data(:,1);
fy = foil.data(:,2);
axis_x = 0.35;
axis_y = 0.034;
angle  = 3.4*pi/180;
Rot = [cos(angle)   sin(angle); ...
       -sin(angle)  cos(angle)];

coords = Rot*[fx'-axis_x; fy'-axis_y];
fx2 = coords(1,:)+axis_x;
fy2 = coords(2,:)+axis_y;
%plot(fx2,fy2, 'b', 'LineWidth', 2);
%xlim([0.6 1.05])
%view(2)


% Reynolds stress
outpos2 = [0.35 0.5 0.4*gr 0.4];
figure(2)
set(gcf,'Units','Normalized')
set(gcf,'OuterPosition',outpos2);
for i=1:3
  if (i==1)     
    var=uu{pr}(iy,:);
    subplot(1,3,i);
    ylabel('$\overline{uu},\overline{vv},\overline{uv}$')
    xlabel('$\alpha$')
   
  elseif (i==2)
    var=vv{pr}(iy,:);
    subplot(1,3,i);
  elseif (i==3)
    var=uv{pr}(iy,:);
    subplot(1,3,i);
  end  

  plot(alpha,var); hold on
  
  tfine = linspace(time(1),time(end),2000);
  vfine = interp1(time,var,tfine);
  ar = areafcn(tfine,vfine,0);
  
  lag0 = 0;
  [oplag area exitflag output] = fminsearch(@(lag) areafcn(tfine,vfine,lag),lag0);
  lag=oplag;
  
  %lag=-1.6;
  afine = getalpha(tfine,lag);
  plot(afine,vfine,'r.');

  if (i==1)
    disp(['Lag for uu: ', num2str(lag,5)])
  elseif (i==2)
    disp(['Lag for vv: ', num2str(lag,5)])
  elseif (i==3)
    disp(['Lag for uv: ', num2str(lag,5)])
  end

  title(['$\phi=',num2str(lag,5),'$'])

  if (i==1)
    ylabel('$\overline{uu}$')
  elseif (i==2)
    ylabel('$\overline{vv}$')
  elseif (i==3)
    ylabel('$\overline{uv}$')
  end  

  xlabel('$\alpha$')

end

if (ifsave)
  SaveFig(gcf,'rs_pr63_ow.eps','plots',1)
end  


% Make a surface plot
j=0;
for i=1:nprofs
   j=j+1;
   bern(:,j,:) = 0.5*(UU{i}+VV{i}+WW{i}) + P{i};
   vort(:,j,:) = Wz{i};
   X(:,j,:)    = x{i};
   Y(:,j,:)    = y{i};
   YN(:,j,:)   = yn{i};
   varuu(:,j,:)= uu{i};
   Um(:,j,:)   = U{i};
   Vm(:,j,:)   = V{i};
end


b1=bern(:,:,ts);
vort1=vort(:,:,ts);
uu1=varuu(:,:,ts);
um1=Um(:,:,ts);
vm1=Vm(:,:,ts);

x1=X(:,:,ts);
y1=Y(:,:,ts);

figure(1)
ss = surf(x1(1:75,:),y1(1:75,:),um1(1:75,:),'EdgeColor','interp');
%contour(x1,y1,abs(vort1), 50, 'EdgeColor','interp')
view(2)
xlim([0.6 1.05])
set(ss,'FaceAlpha', 0.5);
set(ss,'EdgeAlpha', 0.5);
colorbar
if (ifsave)
  SaveFig(gcf,'pts_pr63_ow.eps','plots',1)
end  

%if (ifsave)
%  SaveFig(gcf,'alpha_29_d.png','plots',1)
%end  



ifsave=0;

% outpos3 = [0.10 0.2 0.4*gr 0.7];
% figure(3)
% set(gcf,'Units','Normalized')
% set(gcf,'OuterPosition',outpos3);
% 
% surf(x1,y1,um1,'EdgeColor','interp')
% %contour(x1,y1,abs(vort1), 50, 'EdgeColor','interp')
% view(2)
% xlim([0.4 1.05])
% acycle = getalpha(time(ts),0.);
% title(['$\alpha=',num2str(acycle,3), '$'])
% colorbar
% if (ifsave)
%   SaveFig(gcf,'alpha_29_d.png','plots',1)
% end  
% 
% 
% figure(3)
% surf(x1,y1,uu1,'EdgeColor','interp')
% %contour(x1,y1,abs(vort1), 50, 'EdgeColor','interp')
% view(2)
% xlim([0.4 1.05])
% %tcycle = (time(ts)-ptch_start)/Tosc;
% %title(['Time=',num2str(tcycle,4)])
% 
% acycle = getalpha(time(ts),0.);
% title(['$\alpha=',num2str(acycle,3), '$'])
% colorbar
% 
% if (ifsave)
%   SaveFig(gcf,'alpha_29uu_d.png','plots',1)
% end  




