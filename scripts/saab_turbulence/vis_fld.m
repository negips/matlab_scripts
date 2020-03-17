% Trying to visualize convection

clc 
close all
%warning off all
clear all

%addpath '/home/prabal/workstation/my_matlab_scripts/';

datafiles{1} = 'wing_ff5w5_vO.mat';
datafiles{2} = 'wing_ff5w5_vo.mat';
datafiles{3} = 'wing_ff5w1_vo.mat';
datafiles{4} = 'wing_ff5w2_vo.mat';
datafiles{5} = 'wing_ff5w2_vo_dt2.mat';
datafiles{6} = 'wing_ff2w2_vo.mat';

%datafiles{5} = 'wing_fst.mat';

legs{1} = 'FF=5,W=5 Outflow';
legs{2} = 'FF=5,W=5, dong';
legs{3} = 'FF=5,W=1'; 
legs{4} = 'FF=5,W=2';
legs{5} = 'FF=5,W=2; DT2';
legs{6} = 'FF=2,W=2';


%datafiles{3} = 'les_dxdydz180609_t0225.mat';
%datafiles{4} = 'data_les_kai_1000_k_2_t3.mat';

col = ['kbr'];
mkr = ['o-s'];
fs = 16;


destn='plots/';
pcol=1;
Re = 100000;
nu = 1/Re;
ifnorm = 0;
ifdiff = 0;
iffoil = 0;

xmin=0.15;
xmax=0.70;

ihc(1)=0;
ihc(2)=0;

for ifile=[1 6]

  ifdns=0;    

  les=load(datafiles{ifile});          

  l1=length(les.top);     
  top_positions = zeros(l1,1);
  for i=1:l1
    top_positions(i) = les.top(i).xa;
  end

  l2=length(les.bottom);     
  bot_positions = zeros(l1,1);
  for i=1:l2
    top_positions(i) = les.top(i).xa;
    bot_positions(i) = les.bottom(i).xa;
  end

  for isf=1:1

    h1(ifile,isf)=figure;
    colormap('jet')
%    colormap(flipud(colormap))  

    if (isf==1)
      posi=top_positions;
      lesdata=les.top;
      ttl = 'Top';
    else
      posi=bottom_positions;
      lesdata=les.bottom;
      ttl = 'Bottom';
    end

    l3=length(lesdata);  

    fld      = [];
    surf_x   = [];
    surf_y   = [];  
       
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
%        if ifile~=1  
          lesrt       =      (lesdata(ii).DRTx + lesdata(ii).DRTy + lesdata(ii).DRTz);
%        end

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
        end

        res=lesprod+lesdiss+lesttrns+lesvdiff+lesptrns+lesconv;
        if ~ifdns    
          res = res + lesrt;
        end
        
        dpdyn=gradient(lesdata(ii).P,lesdata(ii).yn);    

        fld1    = lesconv;%dpdyn; %lesconv;

        if (iffoil)
          surf_x1 = lesdata(ii).x;
          surf_y1 = lesdata(ii).y; 
        else
          surf_x1 = posi(ii)*ones(length(fld1),1);
          surf_y1 = lesdata(ii).yp;
        end 
            
        fld    = [fld fld1];
        surf_x = [surf_x surf_x1];
        surf_y = [surf_y surf_y1];    

      end     % if posx

    end       % ii

    if ifdns==1
       dnsfld=fld;
       ttl='Convection';  
    elseif (ifdiff)
       fld=fld-dnsfld;
       ttl='Convection err';
    else
       ttl='Convection';  
    end  
 

    figure(h1(ifile,isf))
    surf(surf_x,surf_y,fld,'EdgeColor','none') 
%    legend(legs,'Interpreter', 'Latex', 'FontSize',fs,'Location', 'Best')
    ylabel('$y_{p}$','Interpreter','Latex','FontSize',fs)
    xlabel('x/C','Interpreter','Latex','FontSize',fs)
    view([0 90])
    ylim([0 500])
    xlim([xmin xmax])  
    colorbar     
    set(gca,'Yscale', 'log')
    title(ttl, 'FontSize', fs)   

    hsave = h1(ifile,isf);
    filename=['cmap_conv' num2str(ifile)];
        
%    SaveFig(hsave, filename, destn, pcol)

  end       % isf


end         % ifile



