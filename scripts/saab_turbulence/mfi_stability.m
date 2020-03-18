clc 
close all
%warning off all
clear all

%addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/';
addpath '/scratch/negi/git_repos/matlabscripts/scripts/';
addpath '/scratch/negi/stability_new/';
%addpath '/scratch/negi/stability_ardeshir/2.9-2.10_PPF/';

%...generate Chebyshev differentiation matrices
%global D0 D1 D2 D3 D4

%N=64; % interpolate onto N+1 chebishev points
%[D0,D1,D2,D3,D4] = ChebMat2(N);

col = ['kbrm'];
mkr = ['--sd'];
fs = 16;
lafs=24;
lgfs=8;
budgetfs = 12;

destn='plots/';
pcol=1;
ifnorm = 0;
ifrtf  = 0;
stab_x = 0.2;


  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

fileind = [2000:2001];
%fileind = [400 410 420 425];
nfiles=length(fileind);
cols_mfi = lines(nfiles);

for mfi= fileind            % 76 first absolute instability
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1; 
  xrange_min = 0.25;
  xrange_max = 0.25;
  
  datafile = sprintf('%s%4.4d%s','re750k_rms/saab750k_rms',mfi,'.mat');
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

  x_eig = 1:20:npos;    
  ecols = jet(length(x_eig));

  ind_pick = [1:100];

  ubtmp = 100;
  uetmp = 1;
  for ix = 1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    dp = gradient(total_p,lesdata(ix).yn);

    surf_x = [surf_x lesdata(ix).x];
    surf_y = [surf_y lesdata(ix).y];
    surf_yn = [surf_y lesdata(ix).yn];
    surf_v = [surf_v lesdata(ix).U];

%----------------------------------------

%   Plot a profile

    indx=find(x_eig==ix);
    if isempty(indx) 
       continue
    end

    [val ind] = min(abs(x_pos - stab_x));
    u_prof = lesdata(ix).U;
    yn_prof = lesdata(ix).yn;
  
    legx{indx}=['$x=',num2str(lesdata(ix).xa,2), '$'];

    figure(4)
    plot(u_prof,yn_prof,'Color', ecols(indx,:)); hold on
    xlabel('$U$')
    ylabel('$y_{n}$')
    ylim([0 0.1])
    legend(legx,'Interpreter','latex', 'FontSize',lgfs, 'Location', 'best')
    grid on
  
    stab_profile = u_prof - u_prof(1);  
  
    yinf = yn_prof(end);
  
    beta=000;
    alpha=100.0;
    alpha_min=10;
    alpha_max=1500;
  
    ifmap=1;
    if (ifmap)
      nalpha=200;
      alpha_range = linspace(alpha_min,alpha_max,nalpha);
    else
      alpha_range= [alpha];
      nalpha=1;
    end
    
    ifOS = 1;
    N=151;
    wr = zeros(N+1,nalpha);
    wi = zeros(N+1,nalpha);
    ar = zeros(N+1,nalpha);
  
    col_os = lines(nalpha);
    uns_ar = [];
    uns_wr = [];
    uns_wi = [];
    mfibreak = 0;
    for ios = 1:nalpha
      alpha_ios = alpha_range(ios);
      [V Eval] = StabilityProperties(stab_profile,yn_prof,N,Re,alpha_ios,beta,ifOS);
      lambdar = real(Eval);
      lambdai = imag(Eval);
  %    figure(300)
  %    plot(lambdar,lambdai, '.', 'Color', col_os(ios,:)); hold on
  %    plot(lambdar/alpha_ios,lambdai, '.', 'Color', col_os(ios,:)); hold on
  %    xlabel('\omega_{r}')
  %    ylabel('\omega_{i}')
  %    grid on
  
      wr(:,ios) = lambdar;
      wi(:,ios) = lambdai;
      ar(:,ios) = alpha_ios*ones(N+1,1);
  
      eig_ind = find(lambdai>0);
      if length(eig_ind)>1
        disp('Multiple unstable eigen-values.')
        mfibreak = 1;
        figure(300)
        plot(lambdar,lambdai, '.', 'Color', col_os(ios,:)); hold on
        plot(lambdar/alpha_ios,lambdai, '.', 'Color', col_os(ios,:)); hold on
        xlabel('$\omega_{r}$')
        ylabel('$\omega_{i}$')
        grid on
       
  %      break
      end  
        
      if ~(isempty(eig_ind))
        uns_wr = [uns_wr lambdar(eig_ind)];
        uns_wi = [uns_wi lambdai(eig_ind)];
        uns_ar = [uns_ar alpha_ios];
      end 
  
    end % ios
  
    if mfibreak
      break
    end  
  
    if length(uns_ar)>0 
      figure(301)
      plot(uns_ar,uns_wi, '-', 'Color', ecols(indx,:)); hold on
      ylabel('$\lambda_{i}$')
      xlabel('$\alpha_{r}$')
      grid on
      legend(legx,'Interpreter','latex', 'FontSize',lgfs, 'Location', 'best')     
 
      figure(7)
      plot(uns_wr,uns_wi, '-', 'Color', ecols(indx,:)); hold on
      xlabel('$\lambda_{r}$', 'FontSize', lafs);
      ylabel('$\lambda_{i}$', 'FontSize', lafs);
      legend(legx,'Interpreter','latex', 'FontSize',lgfs, 'Location', 'best')     

      mfi_lambdar{cnt} = uns_wr;
      mfi_lambdai{cnt} = uns_wi;
      mfi_alphar{cnt}  = uns_ar;
  
    end  

  end       % npos
end   % mfi


figure(1)
colormap jet
surf(surf_x,surf_y,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2)
title(['Time= ' num2str(les.timee)])
% superpose boundary layer height on surface plot
hold on
%plot3(blx_wz,bly_wz,10000*ones(size(blx_wz)), 'k', 'LineWidth', 2)
%view(2)


%% Plot profiles






