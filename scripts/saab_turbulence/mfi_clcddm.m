clc 
close all
%warning off all
clear all

%addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/';
addpath '/scratch/negi/git_repos/matlabscripts/scripts/';
addpath './ardeshir_bl/Evaluate_airfoil'

col = ['kbrm'];
mkr = ['--sd'];
fs = 16;
lafs=24;
budgetfs = 12;

destn='plots/';
pcol=1;
ifnorm = 0;
ifrtf  = 0;
ifprofile=1;
  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

sep_x  = [];
sep_y  = [];
sep_t  = [];

U0=1.0;
chord=1.0;
semichord=chord/2;
k=0.4;
Omega=k*U0/semichord;
Tosc=2*pi/Omega;
ptch_start=6.0;
alpha0 = 3.4;
dalpha=1.0;
phase_shift=-pi/2;
rho=1;
Area=1;
clnorm = 0.5*rho*(U0^2)*Area;

for mfi= [1:1:1025]
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1;    
  
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
     
  top_xa = zeros(l1,1);
  top_ya = zeros(l1,1);
 
  bot_xa = zeros(l2,1);
  bot_ya = zeros(l2,1);
 
  for i=1:l1
    top_xa(i) = les.top(i).xa;
    top_ya(i) = les.top(i).ya;
  end
  
  [s_txa tsort] = sort(top_xa);
  s_tya = top_ya(tsort);
  
  for i=1:l2
    bot_xa(i) = les.bottom(i).xa;
    bot_ya(i) = les.bottom(i).ya;
  end

  [s_bxa bsort] = sort(bot_xa);
  s_bya = bot_ya(bsort);

  lesdata = [les.bottom(bsort(end:-1:1)) les.top(tsort)];

  xall = [s_bxa(end:-1:1); s_txa];
  yall = [s_bya(end:-1:1); s_tya];
  xdiff = [0; diff(xall)];    
  ydiff = [0; diff(yall)];
  dsall = sqrt(xdiff.^2 + ydiff.^2);
  sall = cumsum(dsall);
  
%  figure(101)
%  plot(xall,yall)

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
  surf_P = [];

  pt_x   = [];      
  pt_y   = [];      
  pt_yn  = [];
  pt_v   = []; 

  ww_x   = [];      
  ww_y   = [];      
  ww_yn  = [];
  ww_v   = [];
  wall_p = [];
  wall_sny = [];

  wall_p = zeros(npos,1);
  wall_x = zeros(npos,1);
  wall_y = zeros(npos,1);
  wall_tauw = zeros(npos,1);

  ifsep = 0;

  for ix = 1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    dp = gradient(total_p,lesdata(ix).yn);

    if mfi==1025
      surf_x = [surf_x lesdata(ix).x];
      surf_y = [surf_y lesdata(ix).y];
      surf_yn = [surf_y lesdata(ix).yn];
      surf_v = [surf_v lesdata(ix).uv];
      surf_U = [surf_U lesdata(ix).U];
      surf_P = [surf_P lesdata(ix).P];

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
    end  

    wall_p(ix) = lesdata(ix).P(1);
    wall_sny(ix) = -lesdata(ix).sny(1);   % stored direction in inward normal.
    wall_x(ix) = lesdata(ix).x(1);
    wall_y(ix) = lesdata(ix).y(1);
    wall_tauw(ix) = lesdata(ix).tauw;

  end     % npos
%----------------------------------------

  dx = diff(xall);
  dx = abs([dx(1); dx]);

  p_n_dx = sum(wall_p.*dx.*(sign(wall_sny)'));
  cl = trapz(sall,wall_p.*wall_sny');
  cl_all(mfi) = cl;
  cl_all2(mfi) = p_n_dx;

  cl_time(mfi) = les.timee;
  alpha_all = alpha0 + dalpha*sin(Omega*(cl_time - ptch_start) + phase_shift);

  figure(10)
  plot(alpha_all,cl_all/clnorm); hold on
  plot(alpha_all,cl_all2/clnorm, 'r'); hold on


%  legs{cnt} = ['$t/T_{osc}=' num2str((les.timee-ptch_start)/Tosc) '$'];


end   % mfi


figure(1)
surf(surf_x,surf_y,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2)
title(['$t/T_{osc}=' num2str((les.timee-ptch_start)/Tosc) '$'])








