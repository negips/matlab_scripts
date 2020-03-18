clc 
close all
%warning off all
clear all

%addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/';
%addpath '/scratch/negi/git_repos/matlabscripts/scripts/';

col = ['kbrm'];
mkr = ['--sd'];
fs = 16;
budgetfs = 12;

destn='plots/';
pcol=1;
ifnorm = 0;
ifdstary  = 0;
xin_dstar = 0.5;
iffoil = 1;
  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

for mfi= [1]  %[1:315]

  cnt=cnt+1;    
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

  xrange_min = min(top_xa); % 0.01;
  xrange_max = 0.99;

  isf = 1;

  if (isf==1)  
    ii1 = top_xa>=xrange_min;
    ii2 = top_xa<=xrange_max;
    ii  = find(ii1.*ii2);
    top_trunc = top_xa(ii);
    [x_vals isort] = sort(top_trunc);
  else
    ii1 = bot_xa>=xrange_min;
    ii2 = bot_xa<=xrange_max;
    ii  = find(ii1.*ii2); 
    bot_trunc = bot_xa(ii);
    [x_vals isort] = sort(bot_trunc);
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
  surf_U = [];
  surf_Uinv = [];

  [val dstar_ind] = min(abs(x_vals - xin_dstar));

  ind_pick = [1:100];
  bl_criteria = 0.99;
  
  for ix = 1:npos

%    lesconv     = -1/2*(lesdata(ix).Cxx  + lesdata(ix).Cyy  + lesdata(ix).Czz);
    wz = lesdata(ix).dVdx - lesdata(ix).dUdy;
    wz2 = abs(wz);
    wz_criterion = 0.002*max(wz2);
    wz_ind = find(wz2>wz_criterion);
    wz_ind2 = wz_ind(end);
    wz3 = wz2(wz_ind2:end);
    yn = interp1(wz3,lesdata(ix).yn(wz_ind2:end),wz_criterion);
    delta99_wz(ix) = yn;


%   delta99 using total pressure  
    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    deficit = 1-max(total_p);
    total_p2 = total_p + deficit;

%   boundary layer height            
    yn = interp1(total_p2(ind_pick),lesdata(ix).yn(ind_pick),bl_criteria);
    delta99(ix)=yn;

    yn = interp1(total_p2(ind_pick),lesdata(ix).y(ind_pick),bl_criteria);
    bly(ix)=yn;
    yn = interp1(total_p2(ind_pick),lesdata(ix).x(ind_pick),bl_criteria);
    blx(ix)=yn;

%   free-stream velocity at BL edge  
    Ue(ix) = interp1(total_p2,lesdata(ix).U,bl_criteria); 

%   Beta parameter. Calculated in 2 different ways
    dpdx(ix) = interp1(total_p2,lesdata(ix).dPdx, bl_criteria);
    ind2 = find((lesdata(ix).yn-delta99(ix))<0);
%   Way 1: 
    u1 = lesdata(ix).U;
    u1(ind2) = Ue(ix);
    dstar_y = u1-lesdata(ix).U;
    dstar(ix) = trapz(lesdata(ix).yn,dstar_y);

%   Way 2:      
    u2 = lesdata(ix).U;
    ind2 = find((lesdata(ix).yn-delta99(ix))<0);
    grad_pt = max(ind2)+10;
    if grad_pt>length(lesdata(ix).dUdy)
      grad_pt = max(ind2);
    end  
    dudy = lesdata(ix).dUdy(grad_pt);   % pick dUdy from a point above the boundary layer. 10 is a bit arbitrary. Need to find a better method.
    u2_0 = lesdata(ix).U(max(ind2)+1);
    y2_0 = lesdata(ix).yn(max(ind2)+1);
    u_inv = u2_0 + dudy*(lesdata(ix).yn - y2_0);
    ind3 = find(lesdata(ix).yn>delta99(ix));
    u_inv(ind3) = lesdata(ix).U(ind3); 

    dstar_y2 = u_inv - u2; 
    dstar2(ix) = trapz(lesdata(ix).yn,dstar_y2);

    if (ifdstary) && ix==dstar_ind
      dstar_xin = dstar_y2;          
      figure(30)
      plot(lesdata(ix).yn,dstar_y2)
      legend('\delta_{*}(y)')
    end  


%   Theta using Way 1 
    theta_y = (u1-lesdata(ix).U).*(lesdata(ix).U./u1);
    theta(ix) = trapz(lesdata(ix).yn,theta_y);
    theta_y2 = (u_inv-lesdata(ix).U).*(lesdata(ix).U./u_inv);
    theta2(ix) = trapz(lesdata(ix).yn,theta_y2);

%   Shape Factor using Way 1 
    H(ix)   = dstar(ix)/theta(ix);
    H2(ix)   = dstar2(ix)/theta2(ix);

%   finally beta      
    beta(ix) = dstar(ix)/lesdata(ix).tauw*dpdx(ix);
    beta2(ix) = dstar(ix)/lesdata(ix).tauw*lesdata(ix).dPdx(1);
    beta3(ix) = dstar2(ix)/lesdata(ix).tauw*dpdx(ix);

    redstar(ix)=Ue(ix)*dstar(ix)/nu;
    retheta(ix)=Ue(ix)*theta(ix)/nu;
    retau(ix)=lesdata(ix).ut*delta99(ix)/nu;

%%   Zagarola Smitts scales
    U_zs(ix)      = Ue(ix)*dstar(ix)/delta99(ix);        % velocity scale
    gamma_zs(ix)  = U_zs(ix)/Ue(ix);             % gamma

    dUedx(ix) = interp1(lesdata(ix).yn,lesdata(ix).dUdx,delta99(ix));
    beta_zs(ix) = -(delta99(ix)*dUedx(ix))/(gamma_zs(ix)*Ue(ix));
    pi_zs(ix) = beta_zs(ix)/gamma_zs(ix);


    surf_x = [surf_x lesdata(ix).x(ind_pick)];
    surf_y = [surf_y lesdata(ix).y(ind_pick)];
    surf_yn = [surf_yn lesdata(ix).yn(ind_pick)];
%    surf_v = [surf_v lesconv(ind_pick)];
    surf_v = [surf_v total_p2(ind_pick)];
    surf_U = [surf_U lesdata(ix).U(ind_pick)];
    surf_Uinv = [surf_Uinv u_inv(ind_pick)];

  end       % npos
%----------------------------------------


end   % mfi


figure(1)
surf(surf_x,surf_y,(surf_U),'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2); hold on
if (iffoil)
  plot(top_xa,top_ya, '.b', 'MarkerSize', 4.5)
  plot(bot_xa,bot_ya, '.b', 'MarkerSize', 4.5)
end  

total_p_max = max(surf_v(:));
deficit=1-total_p_max;
surf_v=surf_v+deficit;
%[bl hh] = contour(surf_x,surf_y,surf_v, [bl_criteria bl_criteria],'b');
%blx=bl(1,2:end);
%bly=bl(2,2:end);
blz=100*ones(size(blx));
%delete(hh)
plot3(blx,bly,blz,'k','LineWidth',2)
title(['Time= ' num2str(les.timee)])

figure(2)
plot(surf_x(1,:),delta99, 'LineWidth', 2)
hold on
%plot(x_vals,delta99_wz, '--b')
%plot(surf_x(1,:),dstar, '-k')
%plot(surf_x(1,:),theta, '-r')
%plot(surf_x(1,:),dstar2, '--k')
%plot(surf_x(1,:),theta2, '--r')

grid on
legend({'\delta_{99}', '\delta_{*}', '\theta', '\delta2_{*}', '\theta2' }, 'FontSize', 16, 'Location', 'Best')

figure(3)
plot(surf_x(1,:),retheta)
hold on
plot(surf_x(1,:),redstar, 'r')
legend({'Re_{\theta}', 'Re_{\delta_{*}}'}, 'FontSize', 16, 'Location', 'Best')

figure(4)
plot(surf_x(1,:),beta)
hold on
plot(surf_x(1,:),beta2, 'r')
plot(surf_x(1,:),beta_zs*50, 'k')
ylim([-1 20])
legend({'\beta', '\beta_{2}', '50\beta_{zs}'}, 'FontSize', 16, 'Location', 'Best')

figure(5)
plot(surf_x(1,:),H)
legend({'H'}, 'FontSize', 16)

figure(6)
plot(surf_x(1,:),Ue/4)
hold on
plot(surf_x(1,:),U_zs, 'r')
legend({'U_{e}/4', 'U_{zs}'}, 'FontSize', 16)

%% plot uu_p and bl height at multiple positions



%% plot U and bl height
xin = [];
cols = lines(length(xin));
for ii=1:length(xin)
  [val ind] = min(abs(x_vals - xin(ii)));
%  bl_xin = interp1(surf_x(1,:),delta99,xin(ii));
  bl_xin = delta99(ind);  
  figure(20)
  subplot(3,1,1)
  plot(lesdata(ind).yn/bl_xin,lesdata(ind).U/Ue(ind), 'Color', cols(ii,:))
  hold on
  plot([1 1], [min(lesdata(ind).U/Ue(ind)) max(lesdata(ind).U/Ue(ind))], '--k')
  title(['U/U_{e}'], 'FontSize', 16)
  xlim([0 1.1])
 
  subplot(3,1,2)
  unorm = lesdata(ind).U/U_zs(ind);
  plot(lesdata(ind).yn/bl_xin,unorm, 'Color', cols(ii,:))
  hold on
  plot([1 1], [min(unorm) max(unorm)], '--k')
  title(['U/U_{zs}'], 'FontSize', 16)
  xlim([0 1.2])

  subplot(3,1,3)
  unorm = lesdata(ind).U/(lesdata(ind).ut*retau(ind));
  plts_U(ii)=plot(lesdata(ind).yn/bl_xin,unorm, 'Color', cols(ii,:));
  hold on
  plot([1 1], [min(unorm) max(unorm)], '--k')
  title(['U/U_{zs2}'], 'FontSize', 16)
  xlim([0 1.2])

  leg_x{ii}=['x=' num2str(xin(ii))];
  legend(plts_U,leg_x, 'FontSize', 12)

%% Plot uu
  figure(21)
  subplot(3,1,1)
  uunorm = lesdata(ind).uu/(Ue(ind)^2);
  plot(lesdata(ind).yn/bl_xin,lesdata(ind).uu/Ue(ind), 'Color', cols(ii,:))
  hold on
  ylims = get(gca,'Ylim');
  plot([1 1], ylims, '--k')
  title(['uu/U_{e}^{2}'], 'FontSize', 16)
  xlim([0 1.2])
  
  subplot(3,1,2)
  unorm = lesdata(ind).uu/(lesdata(ind).ut^2);
  plts_uu(ii) = semilogx(lesdata(ind).yp,unorm, 'Color', cols(ii,:));
  hold on
  ylims = get(gca,'Ylim');
  title(['uu^{+}'], 'FontSize', 16)
  xlim([1 1000])
  legend(plts_uu,leg_x, 'FontSize', 12)

  subplot(3,1,3)
  unorm = lesdata(ind).uu/(U_zs(ind)^2);
  semilogx(lesdata(ind).yp,unorm, 'Color', cols(ii,:));
  hold on
  ylims = get(gca,'Ylim');
  title(['uu/U_{zs}^{2}'], 'FontSize', 16)
  xlim([1 1000])


  figure(22)
  plts_vbyu(ii) = plot(lesdata(ind).yn/bl_xin,lesdata(ind).V./lesdata(ind).U, 'Color', cols(ii,:));
  hold on
  xlim([0 1.2])
  legend(plts_vbyu, leg_x, 'FontSize', 12)


end  

%%
%save('delta99.mat', 'delta99', 'x_vals', 'bly', 'blx')

