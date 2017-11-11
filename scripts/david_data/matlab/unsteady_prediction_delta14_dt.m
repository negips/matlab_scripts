% Prediction of cm/cz response based on lag model

clear
clc
close all

addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

ifsave = 1;
iftsave = 0;
destn = 'plots/';
ifcols= 1;
ifdt = 0;

c=0.5;
nu = 1.568E-5;

base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
lfol = length(fol);

run get_files

fs = 18;
lfs = 10;
lw = 1;     % linewidth

%% issues
remove_files = [];                                                          
bad_fits = [];
nostd_shift = [];
steady_shift = [];

%remove_files = [];                      % steady curve seems shifted
%nostd_shift  = [1];
%steady_shift = [0.02];
bad_fits = [22 28];                       % 16,53,59,60 - alpha out of range. 
                                                          % 20 - Some cases seem phase shifted by 90 deg.                                                     
                                                          % 28 - continuous offset in alpha.
                                                          % 43 - does not match


good_fits = -1*[];   % 55 - data is a bit noisy
                                                  % 22 - has spikes in some data points.
                                                  % 33 - eller case 

if (ifdt)
  delays = importdata('d14_delays.txt');
  ndelays =  length(delays.data);
end
                                                  
nfiles=length(filenames);
col1=lines(nfiles);

case_count=0;

kall2=[];
phiall2=[];
gammaall2=[];
betaall2=[];
intgall2=[];
intgbydalpha2=[];
intgnorm2=[];
alpha_all2=[];
theta_all2=[];
ampall2=[];
toffall2=[];
sofstall2=[];
dphiall2=[];
dbetaall2=[];

for jj= [32:45] %1:nfiles

  clc    
  fname=filenames{jj};
  uoo=U0(jj);
  deltacase=defl(1); 
  hfile = [folder fname];
  Re=uoo*c/nu;

  found=1;
%  if case_count>0
%     found=0;   
%     continue;
%  end   
  if deltacase~=14
    found=0;
    continue
  end
  if uoo~=24  && uoo ~= 30
    found=0;
    continue;
  end    

  if ifturb(jj)
    found=0;
    continue;
  end
  if max(remove_files==jj)
    found=0;
    continue;
  end
  if max(bad_fits==jj)
    found=0;
    continue;
  end
%  if max(good_fits==jj)
%    found=0;
%    continue;
%  end
%  if max(nostd_shift==jj)
%    found=0;
%    continue;
%  end

  delayfound = 0;
  deltaT = 0;
  if (ifdt)
    delayind = strfind(delays.textdata(:,1),fname);
    for i=1:ndelays
      if ~isempty(delayind{i})
        deltaT = delays.data(i);
        delayfound = 1;
        break
      end 
    end
  end  

  if ifdt && ~delayfound 
    found=0;    
    continue
  end

%  if (delayfound) 
%    [segments] = split_segments(hfile, uoo, deltacase,deltaT);
%  else    
    [segments] = split_segments(hfile, uoo, deltacase);
%  end  

  flds = fieldnames(segments);
  has_freq = max(strcmp(flds,'rfreq'));
  if ~has_freq
    found=0;    
    continue
  end

  iseg=1;    
  rms_alpha  = rms(segments(iseg).alpha - mean(segments(iseg).alpha));
  mean_alpha = mean(segments(iseg).alpha);
  amp        = (rms_alpha*sqrt(2))*180/pi;
  mean_alpha = mean_alpha*180/pi;

  if mean_alpha>6.0
    disp(['Mean alpha out of range: ' num2str(mean_alpha) '; File no:', num2str(jj)])    
    found=0;
    continue;
  end

  if mean_alpha<1.0
    disp(['Mean alpha out of range: ' num2str(mean_alpha) '; File no:', num2str(jj)])    
    found=0;
    continue;
  end

  if uoo==30       
    model=load('14_static_models_950k.mat');
  elseif uoo==24  
    model=load('14_static_models_765k.mat');
  else
    found=0;
    continue    
  end  
  
  case_count=case_count+1;
  nsegs = length(segments);
  ncases=nsegs; 
%  indicies = [1:nsegs];
%  ncases = length(indicies);
  if (Re<600000)
    re_leg = '565k';
  elseif (Re>600000 && Re < 900000)
    re_leg = '765k';
  else
    re_leg = '950k';
  end

  if max(nostd_shift==jj)
    index=find(nostd_shift==jj);
    steadyphaseshift=steady_shift(index);
    modelalpha = model.alpha + steady_shift(index);
  else    
    steadyphaseshift=0;
    modelalpha = model.alpha + steadyphaseshift;
  end  
  modelcz = model.cz; 
 
  kall=[];
  phiall=[];
  gammaall=[];
  betaall=[];
  intgall=[];
  intgbydalpha=[];
  intgnorm=[];
  alpha_all=[];
  theta_all=[];
  ampall=[];
  toffall=[];
  sofstall=[];
  omgall=[];

%  legs_f{case_count} = [re_leg '; \alpha= ' num2str(mean_alpha), '; No:' num2str(jj)];
  legs_f{case_count} = [re_leg '; \alpha= ' num2str(mean_alpha)];

  dum=1;
  col2=lines(ncases+dum); 
 
  ksegs=0;
  figure(50)
  clf
  alpha_ofst = 0;
  ofst_found = 0;
  for ii = 1:ncases

    iseg=ii;    
    rms_alpha  = rms(segments(iseg).alpha - mean(segments(iseg).alpha));
    mean_alpha = mean(segments(iseg).alpha);
    amp        = (rms_alpha*sqrt(2))*180/pi;
    mean_alpha = mean_alpha*180/pi;
    
    k = segments(iseg).rfreq;
    omega = 2*k*uoo/c;
    f = omega/(2*pi);
    nek_omega = 2*k;
    nek_timeperiod = 2*pi/nek_omega;

    q_time = segments(iseg).qtime;
    q_alpha = segments(iseg).alpha;
  
    p_time = segments(iseg).ptime;
    p_cz = segments(iseg).Cz;
    p_cm = segments(iseg).Cm;
    
    q_cm = interp1(p_time,p_cm,q_time,'pchip');
    q_cz = interp1(p_time,p_cz,q_time,'pchip');

%%  Get offset in alpha
    if (iseg==1)
      if (k<0.02)
        alpha_ofst=0;
        par0(1)=alpha_ofst;
        
        options=optimset('MaxFunEvals',1000,'MaxIter',10000,'TolX',1e-8,'Tolfun',1e-8);
        [par,fval,exitflag,output] = fminsearch(@(par) alpha_offset(par,q_alpha,q_cz,modelalpha,modelcz), par0,options);
        alpha_ofst=par(1);
        disp(['Alpha offset: ' num2str(alpha_ofst*180/pi) ' deg'])
        modelalpha = modelalpha-alpha_ofst*180/pi;
        ofst_found = 1;
      end
    end  

    if (~ofst_found)
      found = 0;
      if case_count>0 && iseg==1
        legs_f(case_count) = [];
        case_count=case_count-1
      end  
      disp(['Can not find steady shift. Min k=' num2str(k) ])
      continue
    end 

    if k<0.05
      continue
    end
%   if jj==22 && iseg==10
%    disp(['Something is wrong in this case: ' fname '; iseg: ' num2str(iseg)])
%    continue;
%   end  

    disp(['File: ', fname, '; No:', num2str(jj)]) 
    disp(['Reduced Frequency: ', num2str(k), ';  iseg= ', num2str(iseg)])
    disp(['Mean alpha: ', num2str(mean_alpha)])
    disp(['Amplitude: ', num2str(amp)])
    disp(['Time Period (nek): ', num2str(nek_timeperiod)])
    legs{ii} = ['Re=', num2str(Re) '; k=', num2str(k)]; 
    disp(['--------------'])
  
    theta=0;
    ini_aoa=1;
    amp2=1*pi/180;
    par0(1)=theta;
    par0(2)=ini_aoa;
    par0(3)=amp2;
    
    options=optimset('MaxFunEvals',1000,'MaxIter',10000,'TolX',1e-8,'Tolfun',1e-8);
    [par,fval,exitflag,output] = fminsearch(@(par) unsteady_alpha2(par,q_time,q_alpha,k,uoo), par0,options);
    theta=par(1);
    mean_alpha2=par(2);
    amp2=par(3);

    if amp2<0
      amp2=abs(amp2);
      theta=theta+pi;
    end

    alpha_pred = mean_alpha2 + amp2*sin(omega*q_time + theta);
    
%    figure(21)
%    plot(q_time,q_alpha*180/pi); hold on
%    plot(q_time,alpha_pred*180/pi, ' ok')
    
    % cm/cz model

%    par0(4)=0;    % steady curve offset      
    options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-8,'TolFun', 1e-8);

%     [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model3(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

%    if (delayfound)
%      phi=-5*pi/180;
%      intg_const=1;
%      toff=0.;
%      par0(1)=phi;
%      par0(2)=intg_const;
%
%      [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model3(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);
%
%      phi=par(1);
%      intg_const=par(2);
%
%    else      
      phi=-5*pi/180;
      intg_const=1;
      toff=0.;
      par0(1)=phi;
      par0(2)=intg_const;
      par0(3)=toff;

      [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model2(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

%      [par,fval,exitflag,output] = fmincon(@(par) unsteady_force_model2(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,[],[],[],[],[-pi 0 -5],[0 10 5],[],options);

      phi=par(1);
      intg_const=par(2);
      toff=par(3);

%    end

%    if (intg_const<0)
%      intg_const=abs(intg_const);
%      phi = phi - pi;
%      if (phi>0)
%        phi=phi-2*pi;
%      end
%    end  
        

%   [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

%    [par,fval,exitflag,output] = fmincon(@(par) unsteady_force_model(par,q_time,q_cz,model.cz,model.alpha,k,uoo,mean_alpha2,amp2,theta),par0,[],[],[],[],[-pi/2 -100],[pi/2 100],[],options);

%    if iseg==1
%      phi=0;
%      toff=0;
%    end  

   
    if (isnan(fval))
      fval    
      continue
    end
    ksegs=ksegs+1;

%    [[mean_alpha2 amp2 theta phi]*180/pi intg_const]
    
    % omega = 2*k*uoo/c;
    time  = q_time+toff;
 
    alpha_pred2 = mean_alpha2 + amp2*sin(omega*time + theta);
    alpha_lagg  = mean_alpha2+amp2*sin(omega*time + theta + phi);
    
%    inst_OMEGA=omega*amp2*cos(omega*time+theta);
    inst_OMEGA=cos(omega*time+theta);
   
    p_motion = intg_const*(inst_OMEGA);
    
    cz_lagg = interp1(modelalpha,modelcz,alpha_lagg*180/pi,'linear');
    
    cz_pred = p_motion + cz_lagg;

    gamma = omega*(toff-deltaT);
    beta  = omega*(toff);

%    [toff deltaT (toff-deltaT) omega gamma gamma*180/pi]

  
    kall=[kall k];
    phiall = [phiall phi];                                  % phase lag quasi-steady term
    gammaall = [gammaall gamma];                            % modified phase gain (-deltaT) added mass
    betaall = [betaall gamma];                              % phase gain added mass
    intgall = [intgall intg_const];
    intgbydalpha=[intgbydalpha intg_const/amp2];
    intgnorm=[intgnorm intg_const/amp2*uoo/30];             % uoo/30 to keep the factor ~ O(1)
    alpha_all=[alpha_all mean_alpha2];
    theta_all=[theta_all theta];
    ampall = [ampall amp2];
    toffall = [toffall toff];                               % added mass time gain 
    omgall = [omgall omega];  
    %    sofstall = [sofstall steady_offset];  

    % Phase plot
    % interpolate onto qtime
   
    figure(10+ii)
    clf
%    subplot(ncases,1,ii)
%    phase_plot = plot(q_alpha*180/pi,q_cz, ' .', 'Color', col2(ii,:)); hold on
    phase_plot = plot(alpha_pred2*180/pi,q_cz, ' .', 'Color', col2(ii,:)); hold on
    ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
    legend(phase_plot, legs(ii), 'Interpreter', 'tex', 'fontsize', lfs, 'Location', 'NorthWest')
%    plot(q_alpha*180/pi,cz_pred, '--', 'Color', col2(ii,:));
%    plot(alpha_pred2*180/pi,cz_pred, '--', 'Color', col2(ncases+dum,:), 'LineWidth', 2);
    plot(alpha_pred2*180/pi,cz_pred, '--o', 'Color', 'k', 'LineWidth', 1);

    if ksegs~=-1
      plot(modelalpha,modelcz,'--k', 'LineWidth', 2)
      xlim([mean_alpha2*180/pi-1.5 mean_alpha2*180/pi+1.5])
    end

    if ifsave 
      filename=['phase_plot_' num2str(jj) '_' num2str(iseg)];
      filename = [re_leg '_' filename];
      SaveFig(gcf,filename, destn, ifcols)
    end

    figure(50)
    time_plot = plot(q_time,q_cz, '. ', 'Color', col2(ii,:), 'LineWidth', 2, 'MarkerSize', 16); hold on
%    time_pred = plot(q_time,cz_pred, '--', 'Color', col2(ncases+dum,:), 'LineWidth', 2);
    time_pred = plot(q_time,cz_pred, '-', 'Color', 'k', 'LineWidth', 2);
    ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('time', 'Interpreter', 'tex', 'FontSize', fs)

    if iftsave
      xlim([min(q_time) max(q_time)])
      ylim([min(q_cz) max(q_cz)])
      figure(50)
      filename=['time_plot_' num2str(jj) '_' num2str(iseg)];
      filename = [re_leg '_' filename];
      SaveFig(gcf,filename, destn, ifcols)
      xlim auto
      ylim auto 
    end


  end

  if found

    if uoo==30
      mkr='-o';
    else
      mkr='-d';
    end 

    if ifsave
      figure(50)
      filename=['time_plot_' num2str(jj)];
      filename = [re_leg '_' filename];
      SaveFig(gcf,filename, destn, ifcols)
    end

    figure(30)
    plot(kall,phiall*180/pi, mkr, 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel('\phi', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on    

    figure(31)
    plot(kall,abs(intgbydalpha), mkr, 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel('Integration constant', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on
    grid on

    figure(32)
    plot(kall,toffall, mkr, 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel('Time offset', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on

    figure(33)
    plot(kall,gammaall*180/pi, mkr, 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel('\gamma', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on    

%    figure(34)
%    plot(kall,abs(intgnorm), mkr, 'Color', col1(jj,:), 'LineWidth', 2)
%%    ylabel('Normalized Integration constant', 'Interpreter', 'tex', 'FontSize', fs)
%    ylabel('C_{1}', 'Interpreter', 'tex', 'FontSize', fs)
%    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
%    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
%    hold on
%    grid on

    figure(35)
    plot(kall,betaall*180/pi, mkr, 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel('\beta', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on    

%    figure(36)
%    plot(kall,abs((toffall+deltaT)./deltaT ), mkr, 'Color', col1(jj,:), 'LineWidth', 2)
%    ylabel('toff/deltaT', 'Interpreter', 'tex', 'FontSize', fs)
%    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
%    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
%    hold on
%    grid on

%    figure(33)
%    plot(kall,sofstall, '-d', 'Color', col1(jj,:))
%    ylabel('Steady offset', 'Interpreter', 'tex', 'FontSize', fs)
%    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
%    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
%    hold on

%   Phase lag derivative
    dphi_all = gradient(phiall,kall);
    dbeta_all = gradient(betaall,kall);

    figure(36)
    plot(kall,dphi_all*180/pi, mkr, 'Color', col1(jj,:), 'LineWidth', 2); hold on
    plot(kall,dbeta_all*180/pi, 'd-', 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel(['\Delta\phi, \Delta\phi'], 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend({'\Delta\phi', '\Delta\beta'}, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on    

    kall2=[kall2 kall];
    phiall2 = [phiall2 phiall];
    gammaall2 = [gammaall2 gammaall];
    betaall2 = [betaall2 betaall];
    intgall2 = [intgall2 intgall];
    intgbydalpha2=[intgbydalpha2 intgbydalpha];
    intgnorm2 = [intgnorm2 intgnorm];
    alpha_all2=[alpha_all2 alpha_all];
    theta_all2=[theta_all2 theta_all];
    ampall2=[ampall2 ampall];
    toffall2=[toffall2 toffall]; 

    dphiall2 = [dphiall2 dphi_all];
    dbetaall2 = [dbetaall2 dbeta_all];

  end

end

% save(['delta' num2str(deltacase) '_predictions.mat'], 'kall2','phiall2','gammaall2','intgall2','intgbydalpha2','intgnorm2','alpha_all2','theta_all2','ampall2','toffall2')


figure(30)
filename=['phase_lag_k-phi.eps'];
%filename = [re_leg '_' filename];
SaveFig(gcf,filename, destn, ifcols)
%
%
%figure(31)
%filename=['phase_lag_k-intg.eps'];
%%filename = [re_leg '_' filename];
%SaveFig(gcf,filename, destn, ifcols)
%
%figure(34)
%filename=['phase_lag_k-intgnorm.eps'];
%%filename = [re_leg '_' filename];
%SaveFig(gcf,filename, destn, ifcols)



