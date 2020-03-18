clc 
close all
%warning off all
clear all

% addpath '/home/prabal/workstation/git_kth/matlabscripts/scripts/';
addpath '/scratch/negi/git_repos/matlabscripts/scripts'

datafiles{1} = 'saab750k001.mat';

%datafiles{5} = 'wing_fst.mat';

legs{1} = '750, N=5';

fileind = [1];
ifref   = [0];

%%
col = ['kbrm'];
mkr = ['--sd'];
fs = 16;
budgetfs = 12;

query_x = 0.7;

destn='plots/';
pcol=1;
%Re = 100000;
%nu = 1/Re;
ifnorm = 1;
ifrtf  = 0;

display(['Normalized: ' num2str(ifnorm)])


ihc(1)=0;
ihc(1)=0;

nfiles = length(fileind);

for ifi=1:nfiles

  ifile=fileind(ifi);
  les=load(datafiles{ifile});
  Re=les.Rer;
  nu = 1/Re;

  if (ifi==1)    
    display(['Re = ' num2str(Re)])
    display(['nu = ' num2str(nu)])
  end    

  ifdns = ifref(ifi);

%  if (strfind(datafiles{ifile},'_vO.mat'))
%    ifdns=1;
%  else
%    ifdns=0;
%  end          

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
      [val ii] = min(abs(top_positions-query_x));
      disp(['Case: ' datafiles{ifile} '; Top profiles at x = ' num2str(top_positions(ii))])
    else
      [val ii] = min(abs(bot_positions-query_x));
      disp(['Case: ' datafiles{ifile} '; Bottom profiles at x = ' num2str(bot_positions(ii))])
    end
      
    ihc(isf)=ihc(isf)+1;   
      
    if isf==1
      if ifdns==1         % dns  
        lesdata=les.top;
      else
        lesdata=les.top;
      end    
      ttl = 'Top';
    elseif isf==2
      if ifdns==1  
        lesdata=les.bottom;
      else  
        lesdata=les.bottom;
      end  
%      dnsdata=dns.bottom;
%      dnsdata = dns.bottom2;
      ttl = 'Bottom';
    end
  
    lesut = lesdata(ii).ut;
    les_budget_norm = (lesut^4)/nu;
    
    ih = 0; 
  % tau_w
  %---------------------------------------------- 
    ih = ih + 1;
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end    
    l1 = length(lesdata);
    for jj=1:l1
      lesutau(jj) = lesdata(jj).tauw;
      wingpos(jj) = lesdata(jj).xa;
    end
    utau_les(ifi,isf) = plot(wingpos,lesutau,['*' col(ifi)], 'MarkerSize', 4); hold on;
    ylabel('$\tau_{w}$', 'Interpreter', 'Latex', 'FontSize', fs)
    xlabel('$x/C$', 'Interpreter', 'Latex', 'FontSize', fs)
    grid on  
  
    title(['$\tau_{w}$ on ' ttl ' Surface'], 'Interpreter', 'Latex', 'FontSize', fs)
    if ifi==nfiles  
      legend(utau_les(:,1),legs(fileind), 'FontSize', fs, 'Location', 'NorthEast', 'Box', 'off')
    end
    ylim([-0.002 0.01])   
  
    hsave = gcf;
    filename=[ttl '_tauw.eps'];
    SaveFig(hsave, filename, destn, pcol)
  
  % delta99
  %---------------------------------------------- 
%    ih = ih + 1;
%    if ihc(isf)==1  
%      h1(isf,ih) = figure;
%    else
%      figure(h1(isf,ih));
%    end    
%    l1 = length(lesdata);
%    cnt = 0;
%    for jj=1:l1
%      loc = lesdata(jj).xa;
%      if loc>=0.20 && loc <=0.999
%        cnt=cnt+1;
%        lesdelta99(cnt) = lesdata(jj).delta99;
% %       dnsdelta99(cnt) = dnsdata(jj).delta99;
%        wingpos2(cnt) = lesdata(jj).xa;
%      end
%    end
%    d99_les(ifi,isf) = plot(wingpos2,lesdelta99,[mkr(ifi) col(ifi)]); hold on;
% %   utau_dns(isf) = plot(wingpos2,dnsdelta99,'sk'); hold on;
%    ylabel('$\delta_{99}$', 'Interpreter', 'Latex', 'FontSize', 16)
%    xlabel('$x/C$', 'Interpreter', 'Latex', 'FontSize', 16)
%  
%    title(['$\delta_{99}$ on ' ttl ' Surface'], 'Interpreter', 'Latex', 'FontSize', 16)
%    legend(d99_les(:,1),{'LES', 'DNS'}, 'FontSize', 16, 'Location', 'West', 'Box', 'off')
%  
%    hsave = gcf;
%    filename=[ttl '_delta99'];
%    SaveFig(hsave, filename, destn, pcol)
  
    %Mean Up
  %---------------------------------------------- 
    ih = ih + 1;
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end    
    upmean_les(ifi,isf) = semilogx(lesdata(ii).yp,lesdata(ii).Up,[mkr(ifi) col(ifi)]); hold on;
    title(['Mean velocity (normalized) (' ttl ') @ x=' num2str(lesdata(ii).xa)], 'FontSize', fs)
    ylabel('$U^{+}$', 'Interpreter', 'Latex')
    xlabel('$y^{+}$', 'Interpreter', 'Latex')
    if ifi==nfiles  
      legend(upmean_les(:,1),legs(fileind), 'FontSize', fs, 'Location', 'Best', 'Box', 'off')
    end  
  
    hsave = gcf;
    filename=[ttl '_meanup'];
    SaveFig(hsave, filename, destn, pcol)

    %Mean U
  %---------------------------------------------- 
    ih = ih + 1;
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end    
    umean_les(ifi,isf) = semilogx(lesdata(ii).yn,lesdata(ii).U,[mkr(ifi) col(ifi)]); hold on;
    title(['Mean velocity (' ttl ') @ x=' num2str(lesdata(ii).xa)], 'FontSize', fs)

    if (ifnorm)  
      ylabel('$U^{+}$', 'Interpreter', 'Latex')
      xlabel('$y^{+}$', 'Interpreter', 'Latex')
    else
      ylabel('$U$', 'Interpreter', 'Latex')
      xlabel('$y$', 'Interpreter', 'Latex')
    end  
    if ifi==nfiles  
      legend(umean_les(:,1),legs(fileind), 'FontSize', fs, 'Location', 'Best', 'Box', 'off')
    end  
  
    hsave = gcf;
    filename=[ttl '_meanu'];
    SaveFig(hsave, filename, destn, pcol)

  
  
    %RMSs
  %---------------------------------------------- 
    ih=ih+1;%2
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end

    if ifnorm
      lesuu = lesdata(ii).uup;   
      lesvv = lesdata(ii).vvp;   
      lesww = lesdata(ii).wwp;   
      lesuv = lesdata(ii).uvp;
      ypos  = lesdata(ii).yp;       
    else
      lesuu = lesdata(ii).uu;   
      lesvv = lesdata(ii).vv;   
      lesww = lesdata(ii).ww;   
      lesuv = lesdata(ii).uv;
      ypos  = lesdata(ii).yn;       
    end   

    rs_les(ifi,isf,1) = semilogx(ypos,lesuu,'-b'); hold on;
    rs_les(ifi,isf,2) = semilogx(ypos,lesvv,'-r'); hold on;
    rs_les(ifi,isf,3) = semilogx(ypos,lesww,'-k'); hold on;
    rs_les(ifi,isf,4) = semilogx(ypos,lesuv,'-m'); hold on;

    if ifdns
      set(rs_les(ifi,isf,1),'LineStyle', 'o')  
      set(rs_les(ifi,isf,2),'LineStyle', 'o')  
      set(rs_les(ifi,isf,3),'LineStyle', 'o')  
      set(rs_les(ifi,isf,4),'LineStyle', 'o')  
    end   

    if (ifnorm)   
      xlabel('$y^{+}$', 'Interpreter', 'Latex')
    else
      xlabel('$y$', 'Interpreter', 'Latex')
    end   

    if ifnorm    
      xlim([0.6 600])
    else
      xlim([0.0001 0.2])  
    end  

    title([ttl ' Reynolds Stresses'], 'FontSize', fs)
  
    hsave = gcf;
    filename=[ttl '_reynolds_stress'];
    SaveFig(hsave, filename, destn, pcol)
  
    %Budget
  %---------------------------------------------- 
  
    lesprod     =  1/2*(lesdata(ii).Pxx  + lesdata(ii).Pyy  + lesdata(ii).Pzz);
    lesdiss     =  1/2*(lesdata(ii).Dxx  + lesdata(ii).Dyy  + lesdata(ii).Dzz);
    lesttrns    =  1/2*(lesdata(ii).Txx  + lesdata(ii).Tyy  + lesdata(ii).Tzz);
    lesptrns    =  1/2*(lesdata(ii).Pixx + lesdata(ii).Piyy + lesdata(ii).Pizz);
    lesvdiff    =  1/2*(lesdata(ii).VDxx + lesdata(ii).VDyy + lesdata(ii).VDzz);
    lesconv     = -1/2*(lesdata(ii).Cxx  + lesdata(ii).Cyy  + lesdata(ii).Czz);
    ypos        = lesdata(ii).yn;    

    if ~ifdns  
      lesrt       =      (lesdata(ii).DRTx + lesdata(ii).DRTy + lesdata(ii).DRTz);
    end  
  
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
      ypos        = lesdata(ii).yp;
    end
  
    res=lesprod+lesdiss+lesttrns+lesvdiff+lesptrns+lesconv;
    if ~ifdns 
      res = res + lesrt;
    end
  
    ih = ih+1;%3
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end

    prod_les(ifi,isf)  = semilogx(ypos, lesprod,  'Color', 'b'             , 'LineStyle', '-', 'LineWidth', 1); hold on;     % production
    diss_les(ifi,isf)  = semilogx(ypos, lesdiss,  'Color', 'r'             , 'LineStyle', '-', 'LineWidth', 1); hold on;     % resolved dissipation
    ttrns_les(ifi,isf) = semilogx(ypos, lesttrns, 'Color', 'm'             , 'LineStyle', '-', 'LineWidth', 1); hold on;
    vdiff_les(ifi,isf) = semilogx(ypos, lesvdiff, 'Color', 'g'             , 'LineStyle', '-', 'LineWidth', 1); hold on;
    ptrns_les(ifi,isf) = semilogx(ypos, lesptrns, 'Color', 'k'             , 'LineStyle', '-', 'LineWidth', 1); hold on;
    conv_les(ifi,isf)  = semilogx(ypos, lesconv,  'Color', [0.68 0.46 0.]  , 'LineStyle', '-', 'LineWidth', 1); hold on;
  
    if ~ifdns && ifrtf 
      diss_tot(ifi,isf)  = semilogx(ypos,lesdiss+lesrt, 'xr', 'LineWidth', 1); hold on;     % total dissipation
      set(diss_tot(ifi,isf),'LineWidth', 1.5)
    end  

    if ifdns
      set(prod_les(ifi,isf),'LineStyle', ' o')  
      set(diss_les(ifi,isf),'LineStyle', ' o')  
      set(ttrns_les(ifi,isf),'LineStyle', ' o')  
      set(vdiff_les(ifi,isf),'LineStyle', ' o')  
      set(ptrns_les(ifi,isf),'LineStyle', ' o')  
      set(conv_les(ifi,isf),'LineStyle', ' o')
    end   
  
    xlabel('$y^{+}$', 'Interpreter', 'Latex')
    title([ttl ' K-Budget @ x= ' num2str(lesdata(ii).xa)])

    if (~ifdns)
      if (ifrtf)
        legend([prod_les(ifi,isf) diss_les(ifi,isf) ttrns_les(ifi,isf) vdiff_les(ifi,isf) ptrns_les(ifi,isf) conv_les(ifi,isf) diss_tot(ifi,isf)], {'Production', 'Resolved Dissipation', 'Turb Transport', 'Viscous Diff.', 'Pressure transport', 'Convection' 'Diss+RTFu'}, 'Location', 'Southeast','FontSize', budgetfs)
      else
        legend([prod_les(ifi,isf) diss_les(ifi,isf) ttrns_les(ifi,isf) vdiff_les(ifi,isf) ptrns_les(ifi,isf) conv_les(ifi,isf)], {'Production', 'Resolved Dissipation', 'Turb Transport', 'Viscous Diff.', 'Pressure transport', 'Convection'}, 'Location', 'Southeast','FontSize', budgetfs)
      end

    end
    if ifnorm    
      xlim([0.6 600])
    else
      xlim([0.0001 0.2])  
    end  
  
    hsave = gcf;
    filename=[ttl '_budget'];
    SaveFig(hsave, filename, destn, pcol)
  
  
    % residual
  %----------------------------------------------  
    ih = ih+1;
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end    
    les_resid(ifi,isf) = semilogx(lesdata(ii).yp,res, [mkr(ifi) col(ifi)]); hold on;
%    dns_resid(isf) = semilogx(dnsdata(ii).yp,res, '-k'); hold on;
    legend([les_resid(1:ifi,isf)], [legs(fileind(1:ifi))], 'FontSize',fs, 'Location', 'Best')  
    title([ttl ' Residuals'], 'FontSize', fs);
    set(gca, 'FontSize', fs)
  
    hsave = gcf;
    filename=[ttl '_residual'];
    SaveFig(hsave, filename, destn, pcol)
  
    % transport summation
  %----------------------------------------------  
    ih = ih+1;
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end    
    les_transport = lesttrns + lesptrns + lesvdiff;
  
    les_tr(ifi) = semilogx(lesdata(ii).yp,les_transport, [mkr(ifi) col(ifi)]); hold on;
    title([ttl ' Net transport'], 'FontSize', fs);
    set(gca, 'FontSize', fs)
  
    hsave = gcf;
    filename=[ttl '_transport'];
    SaveFig(hsave, filename, destn, pcol)


    % Vertical pressure gradient
  %----------------------------------------------  
    ih = ih+1;
    if ihc(isf)==1  
      h1(isf,ih) = figure;
    else
      figure(h1(isf,ih));
    end    

    dpdyn = gradient(lesdata(ii).P,lesdata(ii).yn);
  
    les_tr(ifi) = semilogx(lesdata(ii).yp,dpdyn, [mkr(ifi) col(ifi)]); hold on;
%    title([ttl ' Vertical pressure gradient'], 'FontSize', fs);
    set(gca, 'FontSize', fs)
  
    hsave = gcf;
    filename=[ttl '_dpdy'];
    SaveFig(hsave, filename, destn, pcol)
  
  end   % isf

end   % ifi


