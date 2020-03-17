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
ifprofile=0;
ifplot=0;

kred=0.4;
U0=1.0;
chord=1;
semichord=chord/2;
omega=kred*U0/semichord;
Tosc=2*pi/omega;
ptch_start=6.0;
ptch_amp=1.0;
ini_aoa=3.4;
ini_phase=-pi/2;

  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

uu = [];
vv = [];
ww = [];
uv = [];
xa = [];
ya = [];
yn = [];
T  = [];

x0 = [0.4 0.5 0.6 0.7 0.75 0.8 0.85 0.9];

ic=0;
for mfi= [1:2] % 1560]

   ic=ic+1;
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1;    
  xrange_min = 0.01;
  xrange_max = 1.00;
  
  datafile = sprintf('%s%4.4d%s','re750k_mat/saab750k_pitch',mfi,'.mat');
  les=load(datafile);
  Re=les.Rer;
  nu = 1/Re;
  
  if (mfi==1)    
    display(['Re = ' num2str(Re)])
    display(['nu = ' num2str(nu)])
  end    
  
  l1=length(les.top);
  l2=length(les.bottom);
     
  top_positions = zeros(l1,1);
  bot_positions = zeros(l2,1);
  
  for i=1:l1
    top_positions(i) = les.top(i).xa;
  end
  
  for i=1:l2
    bot_positions(i) = les.bottom(i).xa;
  end
  
  isf = 1;

  if (isf==1)  
    ii1 = top_positions>=xrange_min;
    ii2 = top_positions<=xrange_max;
    ii  = find(ii1.*ii2);
    top_trunc = top_positions(ii);
    [top_trunc isort] = sort(top_trunc);
    allpos = top_trunc;
  else
    ii1 = bot_positions>=xrange_min;
    ii2 = bot_positions<=xrange_max;
    ii  = find(ii1.*ii2); 
    bot_trunc = bot_positions(ii);
    [bot_trunc isort] = sort(bot_trunc);
    allpos = bot_trunc;
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

  surf_x  = [];
  surf_y  = [];
  surf_yn = [];
  surf_yn = [];
  surf_uu = [];
  surf_vv = [];
  surf_ww = [];
  surf_uv = [];

  surf_u = [];
  surf_v = [];

  t=les.timee;
  theta = ini_aoa + ptch_amp*sin(omega*(t-ptch_start)+ini_phase);

  omega=kred*U0/semichord;
  Tosc=2*pi/omega;
  ptch_start=6.0;
  ptch_amp=1.0;
  ini_aoa=3.4;

  indpos = [];
  if ic==1
%   Find points closest to x0s (initial positions)
    for k=1:length(x0)
      [val ind] = min(abs(allpos - x0(k)));
      xa0(k) = lesdata(ind).xa;
      ya0(k) = lesdata(ind).ya;
      indpos = [indpos ind];
    end

    theta0 = ini_aoa + ptch_amp*sin(omega*(t-ptch_start)+ini_phase);
    xa_t = xa0;
    ya_t = ya0;
  else
%   New positions of these points          
    dtheta = (theta-theta0)*pi/180;
    Rot    = [cos(dtheta)   sin(dtheta); ...
              -sin(dtheta)  cos(dtheta)];
    coords = Rot*[xa0; ya0];
    xa_t   = coords(1,:);
    ya_t   = coords(2,:);

%   Nearest index of those positions 
    for k=1:length(xa_t)
      [val ind] = min(abs(allpos - xa_t(k)));
      indpos = [indpos ind];
    end

  end

  for ix = indpos %1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    dp = gradient(total_p,lesdata(ix).yn);

    surf_x =  [surf_x lesdata(ix).x];
    surf_y =  [surf_y lesdata(ix).y];
    surf_yn = [surf_yn lesdata(ix).yn];
    surf_uu = [surf_uu lesdata(ix).uu];
    surf_vv = [surf_vv lesdata(ix).vv];
    surf_ww = [surf_ww lesdata(ix).ww];
    surf_uv = [surf_uv lesdata(ix).uv];

    surf_u = [surf_u lesdata(ix).U];
    surf_v = [surf_v lesdata(ix).V];

  end

  uu{ic} = surf_uu;
  vv{ic} = surf_vv;
  ww{ic} = surf_ww;
  uv{ic} = surf_uv;
  xa{ic} = surf_x;
  ya{ic} = surf_y;
  yn{ic} = surf_yn;

  U{ic}  = surf_u;
  V{ic}  = surf_v;
  T(ic)  = les.timee;
  

  if ifplot
    gray = [0.8 0.8 0.8];
    figure(2)
    plot(surf_uu(1:50,:),surf_yn(1:50,:), 'Color', gray); hold on
    pause(0.01)
  end    


end   % mfi


figure(1)
surf(surf_x,surf_y,surf_uv,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2)
title(['Time= ' num2str(les.timee/Tosc)])

%figure(2)
%plot(tr_time,trx_uv, 'b'); hold on
%plot(tr_time,trx_ww, 'r')
%xlabel('time', 'FontSize', 16)
%ylabel('x/c', 'FontSize', 16)
%title('Transition Location', 'FontSize', 14)
%legend({'u''v''', 'w''w'''}, 'FontSize', 14, 'Location', 'Best')
%SaveFig(gcf,'transition_time.eps', 'plots/', 1)


% save('rmsu.mat', '-v7.3', 'U', 'T', 'Tosc')
% save('rmsv.mat', '-v7.3', 'V', 'T', 'Tosc')
% save('rmsuu.mat', '-v7.3', 'uu','T', 'Tosc')
% save('rmsvv.mat', '-v7.3', 'vv','T', 'Tosc')
% save('rmsww.mat', '-v7.3', 'ww','T', 'Tosc')
% save('rmsuv.mat', '-v7.3', 'uv','T', 'Tosc')
% save('rmsxa.mat', '-v7.3', 'xa','T', 'Tosc')
% save('rmsya.mat', '-v7.3', 'ya','T', 'Tosc')
% save('rmsyn.mat', '-v7.3', 'yn','T', 'Tosc')







