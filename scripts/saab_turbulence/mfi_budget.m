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
  
display(['Normalized: ' num2str(ifnorm)])

cnt=0;

for mfi=1:315
 
  display(['File No: ' num2str(mfi)])

  cnt=cnt+1;    
  xrange_min = 0.01;
  xrange_max = 0.35;
  
  datafile = sprintf('%s%3.3d%s','re100k_mat/re100k_pitch',mfi,'.mat');
  les=load(datafile);
  Re=les.Rer;
  nu = 1/Re;
  
  if (mfi==1)    
    display(['Re = ' num2str(Re)])
    display(['nu = ' num2str(nu)])
  end    
  
  ifdns = 0;
  
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
  
  for isf = 1:1

    if (isf==1)  
      ii1 = top_positions>=xrange_min;
      ii2 = top_positions<=xrange_max;
      ii  = find(ii1.*ii2); 
    else
      ii1 = bot_positions>=xrange_min;
      ii2 = bot_positions<=xrange_max;
      ii  = find(ii1.*ii2); 
    end
      
    if isf==1
      lesdata=les.top(ii);
      ttl = 'Top';
    elseif isf==2
      lesdata=les.bottom(ii);
      ttl = 'Bottom';
    end

    npos = length(lesdata);
    min_val = 1e6;
    max_val = -1e-6;
    tauw = [];
    xvals = [];
    for ix = 1:npos

      total_p = 0.5*(lesdata(ix).UU + lesdata(ix).VV + lesdata(ix).WW) + lesdata(ix).P;
      dp = gradient(total_p,lesdata(ix).yn);
          
      [val ind] = min(lesdata(ix).U);
      if val < min_val
        umin(cnt) = val;
        min_val = val;
        ynmin(cnt) = lesdata(ix).yn(ind);
        xmin(cnt) = lesdata(ix).x(ind);
        ymin(cnt) = lesdata(ix).y(ind);

        [val ind] = max(lesdata(ix).U);
        umax(cnt) = lesdata(ix).U(ind);
        ynmax(cnt) = lesdata(ix).yn(ind);
        xmax(cnt) = lesdata(ix).x(ind);
        ymax(cnt) = lesdata(ix).y(ind);
      end

      [val ind] = max(lesdata(ix).U);
      if val > max_val
        umax_a(cnt) = val;
        max_val = val;
        ynmax_a(cnt) = lesdata(ix).yn(ind);
        xmax_a(cnt) = lesdata(ix).x(ind);
        ymax_a(cnt) = lesdata(ix).y(ind);
      end

      % tauw values
      tauw = [tauw lesdata(ix).tauw];
      xvals = [xvals lesdata(ix).xa];
    
    end
%----------------------------------------
    [xbubble ind] = sort(xvals);
    tauw = tauw(ind);
    ind2 = find(tauw<0,1);
    if ~isempty(ind2) && xbubble(ind2)<0.25
      xbubble2 = xbubble(ind2:end);
      tauw2 = tauw(ind2:end);

      [bubble_start(cnt) ind3] = min(xbubble2);

      xbubble3 = xbubble2(ind3:end);
      tauw3 = tauw2(ind3:end);

%      movtauw3 = sgolayfilt(tauw3,2,5);
%      if length (tauw3)>3
%        [tauw3 long] = movavg(tauw3,3,3,1);
%      end  

      ind4 = find(tauw3>0,1);       % first time tauw gets positive again
      bubble_end(cnt) = xbubble3(ind4);
      bubble_length(cnt) = (bubble_end(cnt) - bubble_start(cnt));
    else
      bubble_length(cnt) = 0;
      bubble_start(cnt) = 0;
      bubble_end(cnt) = 0;
    end        


  
  end   % isf

end   % mfi

figure(1)
plot(abs(umin)); hold on
plot(umax, 'r');
plot(0.20*umax, '--r');
legend({'U_{backflow}', 'U_{max}', '0.2U_{max}'})

figure(2)
plot(ynmin)
title('Location of maximum backflow')

figure(3)
plot(xmin,ymin, '*'); hold on
plot(xmax,ymax, '*r'); hold on

figure(4)
plot(bubble_length)

figure(5)
plot(bubble_start, '*b'); hold on
plot(bubble_end, 'sk');

