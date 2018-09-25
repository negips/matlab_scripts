%% Simulation for the pareto distribution shown by Jordan Peterson

clear
clc
close all

ifsave = 0;

np=1.0e+04;  % No of people simulated
nd=100;     % No of distribution points

N0 = 1.0;   % Initial value;
M_total = N0*np;

p = zeros(np,1);    % variables for each person

nt =1.e+4;             % No of transactions

loss_gain = zeros(np,np);

exchange=1.0E-5;

frames=[];
nframes=1;


h1=figure(1);
set(h1,'Units', 'Normalized');
set(h1,'OuterPosition',[0.03 0.20 0.45 0.55]);

%h2=figure(2);
%set(h2,'Units', 'Normalized');
%set(h2,'OuterPosition',[0.50 0.20 0.45 0.55]);

for i=1:nt
  clc
  i    
  loss_gain = rand(np,np);
  loss_gain = sign(loss_gain - loss_gain');
  cum = sum(loss_gain)*exchange;
  p = p + cum';
  money = N0-p;
  money_percent = money/M_total;
  percent_sorted = sort(money_percent); 

  if mod(i,10)==0
    nframes=nframes+1;
    pends = [min(p) max(p)];
    prange = pends(2)-pends(1);
    phist = (p-pends(1))/prange;
    [ns, xs] = hist(money,nd);

    figure(h1)  
    hist(money,nd);  
    xlim([min(xs) max(xs)])  
    if ifsave  
      fname = sprintf('pareto%5d.eps',nframes);
      SaveFig(gcf,fname,'plots/',1) 
    end  

%    figure(h2)  
%    [ns, xs] = hist(money_percent,nd);
%    hist(money_percent,nd);  
%    xlim([min(xs) max(xs)])

%    cum_sum = ns.*xs;  
%    plot(cum_sum)   
%    if ifsave  
%      fname = sprintf('pareto_cum%5d.eps',nframes);
%      SaveFig(gcf,fname,'plots/',1) 
%    end  
    pause(0.1)

  end      
end


