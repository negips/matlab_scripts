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

%uu = [];
%vv = [];
%ww = [];
%uv = [];
%xa = [];
%ya = [];
%yn = [];
%T  = [];

x0 = [0.2 0.25 0.3 0.35 0.38 0.4:0.01:0.98];

ic=0;
for mfi= [1:79] % 1560]

  ic=ic+1;
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1;    
  xrange_min = 0.01;
  xrange_max = 1.00;
  
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

  npos    = length(lesdata);
  min_val =  1e6;
  max_val = -1e6;
%  tauw    = [];
  xvals   = [];

  surf_x  = [];
  surf_y  = [];
  surf_yn = [];
  surf_yn = [];
  surf_uu = [];
  surf_vv = [];
  surf_ww = [];
  surf_uv = [];

  surf_u  = [];
  surf_v  = [];

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

  if (ic==1)
    nll = [];
    for ix=1:length(indpos)
      x{ix}    = nll;
      y{ix}    = nll;

      xa{ix}   = nll;
      ya{ix}   = nll;
      yn{ix}   = nll;
  
      U{ix}    = nll; 
      V{ix}    = nll;
      P{ix}    = nll;

      UU{ix}   = nll;
      VV{ix}   = nll;
      UV{ix}   = nll;
      WW{ix}   = nll;
      PP{ix}   = nll;
      Wz{ix}   = nll;

      uu{ix}   = nll;
      vv{ix}   = nll;
      ww{ix}   = nll;
      uv{ix}   = nll;
      pp{ix}   = nll;

      ut{ix}   = nll;
      tauw{ix} = nll;

      T{ix}    = nll;

    end        
  end

  cnt2 = 0;
  for ix = indpos %1:npos
    cnt2 = cnt2+1;

%    surf_x =  [surf_x lesdata(ix).x];
%    surf_y =  [surf_y lesdata(ix).y];
%    surf_yn = [surf_yn lesdata(ix).yn];
%    surf_uu = [surf_uu lesdata(ix).uu];
%    surf_vv = [surf_vv lesdata(ix).vv];
%    surf_ww = [surf_ww lesdata(ix).ww];
%    surf_uv = [surf_uv lesdata(ix).uv];
%
%    surf_u = [surf_u lesdata(ix).U];
%    surf_v = [surf_v lesdata(ix).V];

    x{cnt2}    = [x{cnt2} lesdata(ix).x];
    y{cnt2}    = [y{cnt2} lesdata(ix).y];
    xa{cnt2}   = [xa{cnt2} lesdata(ix).xa];
    ya{cnt2}   = [ya{cnt2} lesdata(ix).ya];
    yn{cnt2}   = [yn{cnt2} lesdata(ix).yn];
  
    U{cnt2}    = [U{cnt2} lesdata(ix).U];
    V{cnt2}    = [V{cnt2} lesdata(ix).V];
    P{cnt2}    = [P{cnt2} lesdata(ix).P];

    UU{cnt2}   = [UU{cnt2} lesdata(ix).UU];
    VV{cnt2}   = [VV{cnt2} lesdata(ix).VV];
    UV{cnt2}   = [UV{cnt2} lesdata(ix).UV];
    WW{cnt2}   = [WW{cnt2} lesdata(ix).WW];
    PP{cnt2}   = [PP{cnt2} lesdata(ix).PP];
    Wz{cnt2}   = [Wz{cnt2} (lesdata(ix).dVdx-lesdata(ix).dUdy)];

    uu{cnt2}   = [uu{cnt2} lesdata(ix).uu];
    vv{cnt2}   = [vv{cnt2} lesdata(ix).vv];
    ww{cnt2}   = [ww{cnt2} lesdata(ix).ww];
    uv{cnt2}   = [uv{cnt2} lesdata(ix).uv];
    pp{cnt2}   = [pp{cnt2} (lesdata(ix).PP - lesdata(ix).P.*lesdata(ix).P)];

    ut{cnt2}   = [ut{cnt2} lesdata(ix).ut];
    tauw{cnt2} = [tauw{cnt2} lesdata(ix).tauw];

    T{cnt2}    = [T{cnt2} les.timee];

  end

%  uu{ic} = surf_uu;
%  vv{ic} = surf_vv;
%  ww{ic} = surf_ww;
%  uv{ic} = surf_uv;
%  xa{ic} = surf_x;
%  ya{ic} = surf_y;
%  yn{ic} = surf_yn;
%
%  U{ic}  = surf_u;
%  V{ic}  = surf_v;
%  T(ic)  = les.timee;


end   % mfi


time=T{1};

clearvars -except x y xa ya yn U V P UU VV UV WW PP Wz uu vv ww uv pp ut tauw time kred U0 chord semichord Tosc omega ptch_start ptch_amp ini_aoa ini_phase


save('rms_all.mat', '-v7.3')







