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
ifrtf  = 0;
  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

for mfi= [40]  %[1:315]
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1;    
  xrange_min = 0.1;
  xrange_max = 1.0;
 
%  datafile='wing_dns.mat';
  datafile = sprintf('%s%4.4d%s','re750k_mat20/saab750k_pitch',mfi,'.mat');
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
  surf_U = [];
  surf_V = [];
  surf_Wz = [];


  ind_pick = [1:75];
  wztol = -1e+1; 
 
  for ix = 1:npos

%    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
%    deficit = 1-max(total_p);
%    total_p2 = total_p + deficit;
%    yn = interp1(total_p2(ind_pick),lesdata(ix).yn(ind_pick),0.99);
%    blh(ix)=yn;   
%    lesconv     = -1/2*(lesdata(ix).Cxx  + lesdata(ix).Cyy  + lesdata(ix).Czz);
%
%    dp = gradient(total_p,lesdata(ix).yn);
%    dpdy = gradient(lesdata(ix).P,lesdata(ix).yn);
%
%    wz = lesdata(ix).dVdx - lesdata(ix).dUdy;
%
%    wz2 = abs(wz);
%
%%   Vorticity based boundary layer height      
%    ind2 = find(wz<wztol);
%    ind2 = ind2(end);               % get the last index where vorticity is larger than threshold
%    blh_wz(ix) = interp1(wz(ind2:end), lesdata(ix).yn(ind2:end),wztol);
%    blx_wz(ix) = interp1(wz(ind2:end), lesdata(ix).x(ind2:end),wztol);
%    bly_wz(ix) = interp1(wz(ind2:end), lesdata(ix).y(ind2:end),wztol);
%    blxa_wz(ix)= lesdata(ix).xa;
%    blue_wz(ix) = interp1(wz(ind2:end),lesdata(ix).U(ind2:end),wztol);


    surf_x = [surf_x lesdata(ix).x(ind_pick)];
    surf_y = [surf_y lesdata(ix).y(ind_pick)];
    surf_yn = [surf_y lesdata(ix).yn(ind_pick)];
%    surf_v = [surf_v lesconv(ind_pick)];
%    surf_v = [surf_v total_p2(ind_pick)];
    surf_U = [surf_U lesdata(ix).U(ind_pick)];
    surf_V = [surf_V lesdata(ix).V(ind_pick)];

%    surf_Wz = [surf_Wz wz(ind_pick)];


  end       % npos
%----------------------------------------


end   % mfi


figure(1)
surf(surf_x,surf_y,surf_V,'EdgeColor', 'none', 'FaceColor', 'interp'); colorbar
view(2); hold on

% blz_wz = zeros(size(blx_wz)) + max(surf_Wz(:));
% plot3(blx_wz,bly_wz,blz_wz,'k','LineWidth',2)
% title(['Time= ' num2str(les.timee)])
% 
% figure(2)
% plot(surf_x(1,:),blh_wz)

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

