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

k=0.4;
U0=1.0;
chord=1;
semichord=chord/2;
omega=k*U0/semichord;
Tosc=2*pi/omega;
ptch_start=6.0;
ptch_amp=1.0;
ini_aoa=3.4;

  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

sep_x  = [];
sep_y  = [];
sep_t  = [];


for mfi= [2000:3611]
 
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
  else
    ii1 = bot_positions>=xrange_min;
    ii2 = bot_positions<=xrange_max;
    ii  = find(ii1.*ii2); 
    bot_trunc = bot_positions(ii);
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

  pt_x   = [];      
  pt_y   = [];      
  pt_yn  = [];
  pt_v   = []; 

  ww_x   = [];      
  ww_y   = [];      
  ww_yn  = [];
  ww_v   = [];

  ifsep = 0; 
  
  for ix = 1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    dp = gradient(total_p,lesdata(ix).yn);

    surf_x = [surf_x lesdata(ix).x];
    surf_y = [surf_y lesdata(ix).y];
    surf_yn = [surf_y lesdata(ix).yn];
    surf_v = [surf_v lesdata(ix).uv];

    [v_val ind] = min(lesdata(ix).uv);
    pt_v = [pt_v v_val];
    pt_x = [pt_x lesdata(ix).x(ind)];
    pt_y = [pt_y lesdata(ix).y(ind)];
    pt_yn = [pt_yn lesdata(ix).yn(ind)];

    [v_val ind] = max(lesdata(ix).ww);
    ww_v = [ww_v v_val];
    ww_x = [ww_x lesdata(ix).x(ind)];
    ww_y = [ww_y lesdata(ix).y(ind)];
    ww_yn = [ww_yn lesdata(ix).yn(ind)];

%   Find separation point      
    if ~(ifsep)
      wallshear = lesdata(ix).dUdy(1);
      if wallshear<=0
        sepx = lesdata(ix).xa;
        sepy = lesdata(ix).ya;
        ifsep=1;
      end
    end  

  end       % npos
%----------------------------------------

  % using uv criteria
  [val ind] = min(pt_v);
  uvmin(cnt) = val;
%  pt_v2 = pt_v - val*0.1;
  pt_v2 = pt_v - (-0.0003);
  ind2=find(pt_v2<0,1);
  trx_uv(cnt) = pt_x(ind2);
  tr_uv(cnt)  = pt_v(ind2);

  % using ww criteria
  [val ind] = max(ww_v);
  wwmax(cnt) = val;
%  ww_v2 = ww_v - val*0.1;
  ww_v2 = ww_v - (0.001);
  ind2=find(ww_v2>0,1);
  trx_ww(cnt) = ww_x(ind2);
  tr_ww(cnt)  = ww_v(ind2);

  tr_time(cnt) = les.timee;

%  figure(100)
  %semilogy(top_trunc,abs(pt_v), 'b'); hold on
  %semilogy(top_trunc,ww_v, 'r'); hold off
  
%  plot(top_trunc,abs(pt_v), 'b'); hold on
%  plot(trx_uv(cnt), abs(tr_uv(cnt)), 'sb')
%  plot(top_trunc,ww_v, 'r'); 
%  plot(trx_ww(cnt), tr_ww(cnt), 'or')
%  hold off

  alpha = ini_aoa + ptch_amp*sin(omega*(tr_time - ptch_start)-pi/2);

%  plot(trx_ww(cnt),tr_ww(cnt), 'or'); hold off
%  title(['Time = ' num2str(les.timee)]) 

  if (ifsep)
    sep_x = [sep_x sepx];
    sep_y = [sep_y sepy];
  else
    sep_x = [sep_x lesdata(end).xa];
    sep_y = [sep_y lesdata(end).ya];
  end
  sep_t = [sep_t les.timee];


%  figure(101)
%  plot(alpha,trx_ww, 'r'); hold on
%  plot(alpha,trx_uv, 'b'); hold off
%  ylabel('x/c')
%  xlabel('time')
%
%  figure(102)
%  plot(tr_time,wwmax, 'r'); hold on
%  plot(tr_time,abs(uvmin), 'b'); hold off
%  ylabel('max vals')
%  xlabel('time')



end   % mfi


figure(1)
surf(surf_x,surf_y,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2)
title(['Time= ' num2str(les.timee/Tosc)])

figure(2)
plot(tr_time,trx_uv, 'b'); hold on
plot(tr_time,trx_ww, 'r')
xlabel('time', 'FontSize', 16)
ylabel('x/c', 'FontSize', 16)
title('Transition Location', 'FontSize', 14)
legend({'u''v''', 'w''w'''}, 'FontSize', 14, 'Location', 'Best')
SaveFig(gcf,'transition_time.eps', 'plots/', 1)

figure(3)
plot(tr_time,trx_uv, '-*b'); hold on
plot(tr_time,trx_ww, '-sr')
xlabel('x/c', 'FontSize', 16)
ylabel('time', 'FontSize', 16)
legend({'u''v''', 'w''w'''}, 'FontSize', 14, 'Location', 'Best')
title('Transition Location', 'FontSize', 14)
%SaveFig(gcf,'transition_time.eps', 'plots/', 1)

U0=1.;
kred=0.4;
chord=1.0;
semichord=0.5;
omega=kred*U0/semichord;
ptch_start=6.0;
ptch_amp=1.0;
ini_aoa=3.4;
alpha=ini_aoa + ptch_amp*sin(omega*(tr_time - ptch_start)-pi/2);
figure(4)
plot(alpha,trx_uv, '-b'); hold on
plot(alpha,trx_ww, '-r')
xlabel('\alpha^{o}', 'FontSize', 16)
ylabel('x/c', 'FontSize', 16)
legend({'u''v''', 'w''w'''}, 'FontSize', 14, 'Location', 'Best')
title('Transition Location phase plot', 'FontSize', 14)
SaveFig(gcf,'transition_alpha750k.eps', 'plots/', 1)


figure(5)
plot(sep_t,sep_x, '-b'); hold on
xlabel('x/c', 'FontSize', 16)
ylabel('time', 'FontSize', 16)
legend({'Separation point'}, 'FontSize', 14, 'Location', 'Best')
title('Separation Location', 'FontSize', 14)
SaveFig(gcf,'separation_time.eps', 'plots/', 1)

figure(6)
plot(alpha,sep_x, '* b'); hold on
xlabel('\alpha', 'FontSize', 16)
ylabel('x/c', 'FontSize', 16)
legend({'Separation point'}, 'FontSize', 14, 'Location', 'Best')
title('Separation Location', 'FontSize', 14)
SaveFig(gcf,'separation_phase.eps', 'plots/', 1)


save('tr750k_n9.mat', 'tr_time', 'trx_ww', 'trx_uv', 'alpha')







