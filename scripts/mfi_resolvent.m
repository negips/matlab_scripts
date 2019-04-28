clc 
close all
%warning off all
clear all

addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/stability_ardeshir/';
%addpath '/scratch/negi/git_repos/matlabscripts/scripts/';
%addpath '/scratch/negi/stability_new/';
%addpath '/scratch/negi/stability_ardeshir/2.9-2.10_PPF/';

%...generate Chebyshev differentiation matrices
%global D0 D1 D2 D3 D4

%N=64; % interpolate onto N+1 chebishev points
%[D0,D1,D2,D3,D4] = ChebMat2(N);

col = ['kbrm'];
mkr = ['--sd'];
fs = 16;
lafs=24;
lgfs=20;
budgetfs = 12;

destn='plots/';
pcol=1;
ifnorm = 0;
ifrtf  = 0;
ifbackflow_stab = 1;
stab_x = 0.2;
mkrsiz=12;

  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

bubble_ubmax = [];
bubble_ue = [];
bubble_time = [];
bubble_profile = [];
bubble_yn = [];
bubble_time = [];

fileind = [0];
nfiles=length(fileind);
cols_mfi = lines(nfiles);

for mfi= fileind            % 76 first absolute instability
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1; 
  xrange_min = 0.01;
  xrange_max = 0.90;
  
  datafile = sprintf('%s%4.4d%s','saab750k_pitch',mfi,'.mat');
  les=load(datafile);
  Re=les.Rer;
  nu = 1/Re;
  
  if (cnt==1)    
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
    [x_pos isort] = sort(top_trunc);
  else
    ii1 = bot_positions>=xrange_min;
    ii2 = bot_positions<=xrange_max;
    ii  = find(ii1.*ii2); 
    bot_trunc = bot_positions(ii);
    [x_pos isort] = sort(bot_trunc);
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

  bb_profile = [];
  bb_yn = [];

  ind_pick = [1:100];

  ubtmp = 100;
  uetmp = 1;
  for ix = 1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    dp = gradient(total_p,lesdata(ix).yn);

%    wz = lesdata(ix).dVdx - lesdata(ix).dUdy;

    surf_x = [surf_x lesdata(ix).x];
    surf_y = [surf_y lesdata(ix).y];
    surf_yn = [surf_y lesdata(ix).yn];
    surf_v = [surf_v lesdata(ix).U];

  end       % npos
%----------------------------------------

  for ind=[100] %1:50:npos    

%    ind=length(lesdata);
    u_prof = lesdata(ind).U;
    yn_prof = lesdata(ind).yn;
    
    disp([num2str(ind), ' xa: ', num2str(lesdata(ind).xa, 4)])

    figure(1)
    plot(u_prof,yn_prof,'Color', cols_mfi(cnt,:), 'Marker', '.'); hold on
    xlabel('$U$')
    ylabel('$y_{n}$')
    grid on

    stab_profile = u_prof - u_prof(1);  

    yinf = yn_prof(end);

    beta=000;
    alpha=0.5;
    alpha_min=0.5;
    alpha_max=5;

    if_alphamap=0;
    if (if_alphamap)
      nalpha=10;
      alpha_range = linspace(alpha_min,alpha_max,nalpha);
%      alpha_range = 10.^(linspace(alpha_min,alpha_max,nalpha))
    else
      alpha_range= [alpha];
      nalpha=1;
    end
    
    omega=1.0;
    omega_min=1.0;
    omega_max=50.;
   
    if_omegamap=1;
    if (if_omegamap)
      nomega=5;
      omega_range = linspace(omega_min,omega_max,nomega);
    else
      omega_range= [omega];
      nomega=1;
    end

    ifOS = 0;
    N=250;
    N1=N+1;

        
    cols_ios = lines(nalpha*nomega);
    for ios = 1:nalpha
      for iw = 1:nomega
        omega = omega_range(iw);
        alpha = alpha_range(ios);
        x0 = lesdata(ind).xa;
        y0 = lesdata(ind).ya;
        snx = lesdata(ind).snx;
        sny = lesdata(ind).sny;

        [y1, R, A, B] = PitchingResolvent(stab_profile,yn_prof,x0,y0,snx,sny,N,Re,alpha,beta,omega);
        Rvx = R(1:N1);
        Rvy = R(N1+1:2*N1);
        Rvz = R(2*N1+1:3*N1);      
        Rp  = R(3*N1+1:4*N1);

        cc = (ios-1)*nomega + iw;
        figure(2)
        set(gcf,'Units','normalized');
        set(gcf,'OuterPosition', [0.25 0.35 0.4 0.6]);
        r_r(cc) = plot(real(Rvx),y1,'LineStyle','-' , 'Color', cols_ios(cc,:), 'Marker', '.'); hold on
        r_i(cc) = plot(imag(Rvx),y1,'LineStyle','--', 'Color', cols_ios(cc,:), 'Marker', '.'); hold on
%        r_r(cc) = semilogx(abs(Rvx),yn_prof,'LineStyle','-' , 'Color', cols_ios(cc,:), 'Marker', '.'); hold on
%        r_i(cc) = semilogx(abs(Rvy),yn_prof,'LineStyle','--', 'Color', cols_ios(cc,:), 'Marker', '.'); hold on

        leg_ios{cc} = ['$\alpha=', num2str(alpha), '$; ', '$\omega=',num2str(omega), '$'];
      end   % iw  

    end % ios

    xlabel('$Rx$', 'FontSize', lafs)
    ylabel('$y_{n}$', 'FontSize', lafs)
    grid on
    lg = legend(r_r,leg_ios,'Location','best','FontSize',lgfs); 

    pause(0.01)



  end % ind=1:npos  

end   % mfi


figure(100)
surf(surf_x,surf_y,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2)
title(['Time= ' num2str(les.timee)])
% superpose boundary layer height on surface plot
%hold on
%plot3(blx_wz,bly_wz,10000*ones(size(blx_wz)), 'k', 'LineWidth', 2)
view(2)

%fname = ['lambda_bubble.mat'];
%save(fname, 'mfi_lambdar', 'mfi_lambdai', 'mfi_alphar', 'fileind', 'bubble_profile', 'bubble_ue', 'bubble_yn', 'bubble_time')

%% Plot profiles






