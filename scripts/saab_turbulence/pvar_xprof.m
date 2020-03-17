% Trying to visualize convection

clc 
close all
%warning off all
clear all

addpath '/home/prabal/workstation/my_matlab_scripts/';

datafiles{1} = 'wing_ff5w5_vO.mat';
datafiles{2} = 'wing_ff5w5_vo.mat';
datafiles{3} = 'wing_ff5w1_vo.mat';
datafiles{4} = 'wing_ff5w2_vo.mat';
datafiles{5} = 'wing_ff5w2_vo_dt2.mat';
datafiles{6} = 'wing_ff2w2_vo.mat';
datafiles{7} = 'wing_ff1w2_vo.mat';

%datafiles{5} = 'wing_fst.mat';

legs{1} = 'FF=5,W=5 Outflow';
legs{2} = 'FF=5,W=5, dong';
legs{3} = 'FF=5,W=1'; 
legs{4} = 'FF=5,W=2';
legs{5} = 'FF=5,W=2; DT2';
legs{6} = 'FF=2,W=2';
legs{7} = 'FF=1,W=2';

%postfix = {'vO', 'vo' 'vo' 'vo' 'vo' 'vo'};
postfix = legs;

fileind = [1 6];

mkr = ['osdx'];
fs = 16;
ms = 7;

destn='plots/';
pcol=1;
Re = 100000;
nu = 1/Re;
ifnorm = 1;

xmin = 0.15;
xmax = 0.50;

yppts=[12 25 50];
lyp = length(yppts);

COord=lines(lyp);

ihc(1)=0;
ihc(2)=0;


h1=figure;
set(gca,'ColorOrder',COord)

for ifi=1:length(fileind)

  ifdns=0;    
  ifile=fileind(ifi) 

  les=load(datafiles{ifile}); 

  l1=length(les.top);     
  top_positions = zeros(l1,1);
  for i=1:l1
    top_positions(i) = les.top(i).xa;
  end

  l2=length(les.bottom);     
  bot_positions = zeros(l2,1);
  for i=1:l2
    bot_positions(i) = les.bottom(i).xa;
  end
 
  for isf=1:1

    if (isf==1)
      posi=top_positions;
      lesdata=les.top;
      ttl = 'Top';
      l3=l1;
    else
      posi=bot_positions;
      lesdata=les.bottom;
      ttl = 'Bottom';
      l3=l2;
    end

    pvar_val = [];
    pvar_x   = [];  
       
    for ii=1:l3
      
      if posi(ii)>=xmin && posi(ii)<=xmax
  
        lesut = lesdata(ii).ut;
        les_budget_norm = (lesut^4)/nu;

        lesprod     =  1/2*(lesdata(ii).Pxx  + lesdata(ii).Pyy  + lesdata(ii).Pzz);
        lesdiss     =  1/2*(lesdata(ii).Dxx  + lesdata(ii).Dyy  + lesdata(ii).Dzz);
        lesttrns    =  1/2*(lesdata(ii).Txx  + lesdata(ii).Tyy  + lesdata(ii).Tzz);
        lesptrns    =  1/2*(lesdata(ii).Pixx + lesdata(ii).Piyy + lesdata(ii).Pizz);
        lesvdiff    =  1/2*(lesdata(ii).VDxx + lesdata(ii).VDyy + lesdata(ii).VDzz);
        lesconv     = -1/2*(lesdata(ii).Cxx  + lesdata(ii).Cyy  + lesdata(ii).Czz);
        if (~ifdns)  
          lesrt     =      (lesdata(ii).DRTx + lesdata(ii).DRTy + lesdata(ii).DRTz);
        end
        tke         = (lesdata(ii).uu + lesdata(ii).vv + lesdata(ii).ww)/2;

        if (ifnorm)
          lesprod     =   lesprod /les_budget_norm; 
          lesdiss     =   lesdiss /les_budget_norm;
          lesttrns    =   lesttrns/les_budget_norm;
          lesptrns    =   lesptrns/les_budget_norm;
          lesvdiff    =   lesvdiff/les_budget_norm;
          lesconv     =   lesconv /les_budget_norm;
          if ~ifdns    
            lesrt     =   lesrt   /les_budget_norm;
          end
          tke         =   tke     /les_budget_norm;        
        end

        res=lesprod+lesdiss+lesttrns+lesvdiff+lesptrns+lesconv;
        if ~ifdns    
          res = res + lesrt;
        end

        pvar_val1 = [];
        pvar_x1 = [];

%       Assign variable                
        pvar = lesconv;
          
        for kk=1:lyp
          [val jj] = min(abs(lesdata(ii).yp - yppts(kk)));
          pvar_val1 = [pvar_val1 pvar(jj)];
          pvar_x1   = [pvar_x1 posi(ii)];
          legs{ifi,isf,kk} = ['$y^{+}=' num2str(lesdata(ii).yp(jj)) '$ ', postfix{ifile}]; 
        end
            
        pvar_val = [pvar_val; pvar_val1];
        pvar_x = [pvar_x; pvar_x1];

      end     % if posx

    end       % ii 

    figure(h1)
    plot(pvar_x,pvar_val, 'LineStyle',mkr(ifi)); 
    hold on
    legend(legs(ifi,isf,:),'Interpreter', 'Latex', 'FontSize',fs,'Location', 'Best')
    ylabel('Convection', 'FontSize', fs)
    xlabel('x/C', 'FontSize', fs)

    if ifi>0  
      h2(ifi,isf)=figure;
    else
      figure(h2)
    end          
    subplot(lyp,1,1)
%    set(gca,'ColorOrder',COord);    
    for kk=1:lyp
      subplot(lyp,1,kk)
      sc2(ifi,isf,kk)=plot(pvar_x(:,kk),pvar_val(:,kk), 'x','Color', COord(kk,:));
%      set(sc2(ifi,isf,kk),'Color',COord(kk,:))
  %    ylim([-5 5])  
      legend({['$y^{+}=' num2str(yppts(kk)) '$ ' postfix{ifile}]},'Interpreter', 'Latex', 'FontSize',fs,'Location', 'SouthWest','Box','off')
      hold on
    end        

            
    hsave = h2(1,isf);
    filename=['conv_x' num2str(ifi)];
    SaveFig(hsave, filename, destn, pcol)
  
  end       % isf


end         % ifi

hsave = h1;
filename=['conv_x'];
SaveFig(hsave, filename, destn, pcol)





