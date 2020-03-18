clc 
close all
%warning off all
clear all

%addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/';
addpath '/scratch/negi/git_repos/matlabscripts/scripts/';

col = ['kbrm'];
mkr = ['--sd'];
fs = 16;
budgetfs = 12;

destn='plots/';
pcol=1;
ifnorm = 0;
ifrtf  = 0;
iffoil = 1;
  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

for mfi= [1]  %[1:315]
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1;    
  xrange_min = 0.01;
  xrange_max = 0.98;
 
  datafile = sprintf('%s%3.3d%s','saab750k',mfi,'.mat');
%  datafile='wing_dns.mat';
  les=load(datafile);
  Re=les.Rer;
  nu = 1/Re;
 
  display(['File: ' datafile])
 
  if (mfi==1)    
    display(['Re = ' num2str(Re)])
    display(['nu = ' num2str(nu)])
  end    
  
  l1=length(les.top);
  l2=length(les.bottom);
     
  top_xa = zeros(l1,1);
  bot_xa = zeros(l2,1);

  top_ya = zeros(l1,1);
  bot_ya = zeros(l2,1);
 
  for i=1:l1
    top_xa(i) = les.top(i).xa;
    top_ya(i) = les.top(i).ya;
  end
  
  for i=1:l2
    bot_xa(i) = les.bottom(i).xa;
    bot_ya(i) = les.bottom(i).ya;
  end
  
  isf = 1;

  if (isf==1)  
    ii1 = top_xa>=xrange_min;
    ii2 = top_xa<=xrange_max;
    ii  = find(ii1.*ii2);
    top_trunc = top_xa(ii);
    [top_trunc isort] = sort(top_trunc);
  else
    ii1 = bot_xa>=xrange_min;
    ii2 = bot_xa<=xrange_max;
    ii  = find(ii1.*ii2); 
    bot_trunc = bot_xa(ii);
    [bot_trunc isort] = sort(bot_trunc);
  end
    
  if isf==1
    lesdata=les.top(ii);
    lesdata=lesdata(isort);
    ttl = 'Top';
  elseif isf==2
    lesdata=les.bottom(ii);
    lesdata=lesdata(isort);
    ttl = 'Bottom';
  end

  npos = length(lesdata);
  min_val = 1e6;
  max_val = -1e6;
  tauw = [];
  xvals = [];

  surf_x = [];
  surf_y = [];
  surf_yn = [];
  surf_v = [];
  ind_pick = [1:100];

  bl_criteria = 0.99;
  
  for ix = 1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    deficit = 1-max(total_p);
    total_p2 = total_p + deficit;
    yn = interp1(total_p2(ind_pick),lesdata(ix).yn(ind_pick),bl_criteria);
    delta99(ix)=yn;   

    blx(ix) = interp1(lesdata(ix).yn,lesdata(ix).x,yn);  
    bly(ix) = interp1(lesdata(ix).yn,lesdata(ix).y,yn);  

%    lesconv     = -1/2*(lesdata(ix).Cxx  + lesdata(ix).Cyy  + lesdata(ix).Czz);
    dpdx(ix) = interp1(total_p2(ind_pick),lesdata(ix).dPdx, bl_criteria);
    Ue = interp1(total_p2(ind_pick),lesdata(ix).U, bl_criteria);  
    u1 = lesdata(ix).U-Ue;
    ind2 = find(u1<0);
    u1 = lesdata(ix).U;
    u1(ind2) = Ue;
    dstar_y = u1-lesdata(ix).U;
    dstar(ix) = trapz(lesdata(ix).yn,dstar_y);
    theta_y = (u1-lesdata(ix).U).*(lesdata(ix).U./u1);
    theta(ix) = trapz(lesdata(ix).yn,theta_y);
    H(ix)   = dstar(ix)/theta(ix);

    beta(ix) = dstar(ix)/lesdata(ix).tauw*dpdx(ix);
    beta2(ix) = dstar(ix)/lesdata(ix).tauw*lesdata(ix).dPdx(1);

    redstar(ix)=Ue*dstar(ix)/nu;
    retheta(ix)=Ue*theta(ix)/nu;
    retau(ix)=lesdata(ix).ut*delta99(ix)/nu;

    dp = gradient(total_p,lesdata(ix).yn);
    dpdy = gradient(lesdata(ix).P,lesdata(ix).yn);

    surf_x = [surf_x lesdata(ix).x(ind_pick)];
    surf_y = [surf_y lesdata(ix).y(ind_pick)];
    surf_yn = [surf_y lesdata(ix).yn(ind_pick)];
%    surf_v = [surf_v lesconv(ind_pick)];
    surf_v = [surf_v total_p2(ind_pick)];
%    surf_v = [surf_v lesdata(ix).U(ind_pick)];

  end       % npos
%----------------------------------------


end   % mfi


figure(1)
surf(surf_x,surf_y,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2); hold on
if (iffoil)
  plot(top_xa,top_ya, '.')
  plot(bot_xa,bot_ya, '.')
end  

total_p_max = max(surf_v(:));
deficit=1-total_p_max;
surf_v=surf_v+deficit;
%[bl hh] = contour(surf_x,surf_y,surf_v, [bl_criteria bl_criteria],'b');
%blx=bl(1,2:end);
%bly=bl(2,2:end);
blz=1000*ones(size(blx));
%delete(hh)
plot3(blx,bly,blz,'k','LineWidth',2)
title(['Time= ' num2str(les.timee)])

figure(2)
plot(surf_x(1,:),delta99)
hold on
plot(surf_x(1,:),dstar, '-k')
plot(surf_x(1,:),theta, '-r')
grid on
legend({'\delta', '\delta_{*}', '\theta'}, 'FontSize', 16, 'Location', 'Best')

figure(3)
plot(surf_x(1,:),retheta)
hold on
plot(surf_x(1,:),redstar, 'r')
legend({'Re_{\theta}', 'Re_{\delta_{*}}'}, 'FontSize', 16, 'Location', 'Best')

figure(4)
plot(surf_x(1,:),beta)
hold on
%plot(surf_x(1,:),beta2, 'r')
ylim([-20 20])


%% plot U and bl height
xin = 0.7;
[val ind] = min(abs(top_trunc - xin));
bl_xin = interp1(surf_x(1,:),delta99,xin);

figure(5)
plot(lesdata(ind).yn,lesdata(ind).U)
hold on
plot([bl_xin bl_xin], [min(lesdata(ind).U) max(lesdata(ind).U)], '--k')



%figure(1)
%plot(abs(umin)); hold on
%plot(umax, 'r');
%plot(0.15*umax, '--r');
%
%
%figure(2)
%plot(ynmin)
%
%figure(3)
%plot(xmin,ymin, '*'); hold on
%plot(xmax,ymax, '*r'); hold on
%
%figure(4)
%plot(bubble_length)
%
%figure(5)
%plot(bubble_start, 'b'); hold on
%plot(bubble_end, 'k');

