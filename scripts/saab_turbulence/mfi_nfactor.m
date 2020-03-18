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


for mfi= [10 450]
 
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

  wall_p = zeros(npos,1);
  wall_x = zeros(npos,1);
  wall_y = zeros(npos,1);
  wall_tauw = zeros(npos,1);

  ifsep = 0; 
  
  for ix = 1:npos

    total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
    dp = gradient(total_p,lesdata(ix).yn);

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

%   Find separation point      
    if ~(ifsep)
      wallshear = lesdata(ix).dUdy(1);
      if wallshear<=0
        sepx = lesdata(ix).xa;
        sepy = lesdata(ix).ya;
        ifsep=1;
      end
    end  
    wall_p(ix) = lesdata(ix).P(1);
    wall_x(ix) = lesdata(ix).x(1);
    wall_y(ix) = lesdata(ix).y(1);
    wall_tauw(ix) = lesdata(ix).tauw;

  end     % npos
%----------------------------------------

  % using uv criteria
  [val ind] = min(pt_v);
  pt_v2 = pt_v - val*0.1;
  ind2=find(pt_v2<0,1);
  trx_uv(cnt) = pt_x(ind2);
  tr_uv(cnt)  = pt_v(ind2);

  % using ww criteria
  [val ind] = max(ww_v);
  ww_v2 = ww_v - val*0.1;
  ind2=find(ww_v2>0,1);
  trx_ww(cnt) = ww_x(ind2);
  tr_ww(cnt)  = ww_v(ind2);

  tr_time(cnt) = les.timee;

  Pref=1e5;
  Tref=300;
  Rgas=287.15;
  Rhoref = Pref/(Tref*Rgas);
  Re=7.5e5;
  Visc=visc(Tref);
  chord=1;
  Uref=Re*Visc/Rhoref/chord;
  Mach=Uref/sqrt(Rgas*1.4*Tref);

  param.Re = Re;
  param.Mach = Mach;
  param.Tref = Tref;
  param.fmax = 5e-3/(2*pi*Visc/Rhoref/Uref^2)*1; % just a guess
  sweep=0;
  param.sweep = sweep;
  param.chord = chord;
  Ntr=8;

  result = GetNFac(wall_x,wall_y,wall_p,param,1);
  figure(2)
  plot(result.xn,result.nts); hold on
  ylabel('n-factor')

  figure(3)
  plot(wall_x,wall_p); hold on
  ylabel('pressure')

  figure(4)
  plot(wall_x,wall_tauw); hold on
  ylabel('$\tau_{w}$')
 
end   % mfi


figure(1)
surf(surf_x,surf_y,surf_v,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceColor', 'interp'); colorbar
view(2)
title(['Time= ' num2str(les.timee)])








