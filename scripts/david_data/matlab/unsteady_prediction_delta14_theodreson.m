% Prediction of cm/cz response based on lag model

clear
clc
close all

addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

ifsave = 1;
iftsave = 0;
destn = 'plots/';
ifcols= 1;
ifdt = 1;

c=0.5;
nu = 1.568E-5;

base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
lfol = length(fol);

run get_files

fs = 16;
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
intgall2=[];
intgbydalpha2=[];
intgnorm2=[];
alpha_all2=[];
theta_all2=[];
ampall2=[];
toffall2=[];
sofstall2=[];

for jj= [10] %1:nfiles

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

%    k0=0;
%    kmax=1000;
%    nk=kmax*100;
%    k=linspace(k0,kmax,nk);
    theoderson_k = k;
    theoderson_C = besselh(1,2,k)./(besselh(1,2,k) + i*besselh(0,2,k));
    theoderson_F=real(theoderson_C);
    theoderson_G=imag(theoderson_C);
    theoderson_U0=1;
    theoderson_chord=0.5;
    theoderson_b=theoderson_chord/2;
    theoderson_omega_k = k*theoderson_U0/theoderson_b;
    theoderson_a=-0.15/2;
%    theoderson_ptch_amp=1.0*pi/180;
%    theoderson_r2 = sqrt(1 + (theoderson_a*theoderson_omega_k).^2);
%    theoderson_added_m_eff_amp = -pi*theoderson_ptch_amp*theoderson_omega_k.*theoderson_r2;
%    theoderson_added_m_phi = asin(1/theoderson_r2);


    disp(['File: ', fname, '; No:', num2str(jj)]) 
    disp(['Reduced Frequency: ', num2str(k), ';  iseg= ', num2str(iseg)])
    disp(['Mean alpha: ', num2str(mean_alpha)])
    disp(['Amplitude: ', num2str(amp)])
    disp(['Time Period (nek): ', num2str(nek_timeperiod)])
    legs{ii} = ['Re=', num2str(Re) '; k=', num2str(k)]; 
    disp(['--------------'])

    % Get mean angle, pitching amplitude and initial phase
    theta=0;
    ini_aoa=1;
    amp2=1*pi/180;
    par0(1)=theta;
    par0(2)=ini_aoa;
    par0(3)=amp2;
    
    options=optimset('MaxFunEvals',1000,'MaxIter',10000,'TolX',1e-8,'Tolfun',1e-8);
    [par,fval,exitflag,output] = fminsearch(@(par) unsteady_alpha2(par,q_time,q_alpha,k,uoo), par0,options);

    theta=par(1);             % inital phase
    mean_alpha2=par(2);       % Mean alpha
    amp2=par(3);              % Pitching amplitude

    % Keep amplitude positive
    if amp2<0
      amp2=abs(amp2);
      theta=theta+pi;
    end

    % modelled alpha
    alpha_pred = mean_alpha2 + amp2*sin(omega*q_time + theta);
    
    % cm/cz model

%    par0(4)=0;    % steady curve offset      
    options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-8,'TolFun', 1e-8);

%   Reduced mass contribution is now precalculated up to a phase shift.    
    phi=-0*pi/180;
    par0(1)=phi;              % phase shift for both terms

    [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_theoderson_model(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

    phi=par(1);
   
    if (isnan(fval))
      fval    
      continue
    end
    ksegs=ksegs+1;

%    [[mean_alpha2 amp2 theta phi]*180/pi intg_const]
    
    % omega = 2*k*uoo/c;
    time  = q_time;
 
    alpha_pred2 = mean_alpha2 + amp2*sin(omega*time + theta);
    [cz_pred] = rebuild_unsteady_force(time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta,phi);

  
    kall=[kall k];
    phiall = [phiall phi];
%    gammaall = [gammaall gamma];
%    intgall = [intgall intg_const];
%    intgbydalpha=[intgbydalpha intg_const/amp2];
%    intgnorm=[intgnorm intg_const/amp2*uoo/30];             % uoo/30 to keep the factor ~ O(1)
    alpha_all=[alpha_all mean_alpha2];
    theta_all=[theta_all theta];
%    ampall = [ampall amp2];
%    toffall = [toffall toff-deltaT]; 
%    omgall = [omgall omega];  
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
    time_plot = plot(q_time,q_cz, '-', 'Color', col2(ii,:), 'LineWidth', 2); hold on
%    time_pred = plot(q_time,cz_pred, '--', 'Color', col2(ncases+dum,:), 'LineWidth', 2);
    time_pred = plot(q_time,cz_pred, ' o', 'Color', 'k', 'LineWidth', 1);
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

    figure(34)
    plot(kall,abs(intgnorm), mkr, 'Color', col1(jj,:), 'LineWidth', 2)
%    ylabel('Normalized Integration constant', 'Interpreter', 'tex', 'FontSize', fs)
    ylabel('C_{1}', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on
    grid on

    figure(35)
    plot(kall,abs((toffall+deltaT)./deltaT ), mkr, 'Color', col1(jj,:), 'LineWidth', 2)
    ylabel('toff/deltaT', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on
    grid on

%    figure(33)
%    plot(kall,sofstall, '-d', 'Color', col1(jj,:))
%    ylabel('Steady offset', 'Interpreter', 'tex', 'FontSize', fs)
%    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
%    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
%    hold on
    kall2=[kall2 kall];
    phiall2 = [phiall2 phiall];
    gammaall2 = [gammaall2 gammaall];
    intgall2 = [intgall2 intgall];
    intgbydalpha2=[intgbydalpha2 intgbydalpha];
    intgnorm2 = [intgnorm2 intgnorm];
    alpha_all2=[alpha_all2 alpha_all];
    theta_all2=[theta_all2 theta_all];
    ampall2=[ampall2 ampall];
    toffall2=[toffall2 toffall]; 

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
figure(34)
filename=['phase_lag_k-intgnorm.eps'];
%filename = [re_leg '_' filename];
SaveFig(gcf,filename, destn, ifcols)



