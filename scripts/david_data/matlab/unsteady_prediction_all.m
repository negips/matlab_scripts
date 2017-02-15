% Prediction of cm/cz response based on lag model

clear
clc
close all

destn = 'plots/';
ifcols= 1;

%model.alpha=model.alpha;

c=0.5;
nu = 1.568E-5;

base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
lfol = length(fol);

fs = 16;    % fontsize
lw = 1;     % linewidth

[status,result] = system(['ls ' folder '*']);

inds1 = strfind(result,fol);
inds2 = strfind(result,'.h5');

nfiles = length(inds2);
U0=zeros(nfiles,1);
alpha=zeros(nfiles,1)-99;
defl=zeros(nfiles,1)-99;
ifturb=zeros(nfiles,1);

disp(['N files: ', num2str(nfiles)])

for i=1:nfiles

    ind1 = inds1(i)+lfol;
    ind2 = inds2(i)+2;
    fname=result(ind1:ind2);
    filenames{i}=fname;

    inds4=strfind(fname,'alphasweep');
    inds3=strfind(fname,'_');

    if ~isempty(inds4)
      ind1=2;
      ind2=inds3(1)-1;      
      U0(i) = str2double(fname(ind1:ind2));
      ind1=inds3(1)+2;
      ind2=inds3(2)-1;
      defl(i) = str2double(fname(ind1:ind2));
      if isnan(defl(i)) && i>1 
        defl(i)=defl(i-1);
      end  
      continue
    end
   
    if (length(inds3)>2)
      % Freestream
      ind1=2;
      ind2=inds3(1)-1;
      U0(i) = str2double(fname(ind1:ind2));
      % apha
      ind1=inds3(1)+2;
      ind2=inds3(2)-1;
      alpha(i) = str2double(fname(ind1:ind2));
      % flap deflection      
      ind1=inds3(2)+2;
      ind2=inds3(3)-1;
      defl(i) = str2double(fname(ind1:ind2));
    else
      % Freestream
      ind1=2;
      ind2=inds3(1)-1;
      U0(i) = str2double(fname(ind1:ind2));
      ind1=inds3(1)+1;
      ind2=inds3(2)-1;
      % flap deflection
      defl(i) = str2double(fname(ind1:ind2));
    end

    inds4=strfind(fname,'turb');
    if ~isempty(inds4)
      ifturb(i)=1;
    end  
    
end

fs = 16;
lfs = 12;

%% issues
remove_files = [-22 29 31 32 35 37];                      % steady curve seems shifted
nostd_shift  = [10  12 14 22  41   45   47 49   51    52 54    55];
steady_shift = [0.1 0  0  0.1 0.08 -0.8 0  1.8  -0.75 0  -0.8  0.25];
bad_fits = [16 20 43 53 59 60];                       % 16,53,59,60 - NaNs, 20 - something funny here.
good_fits = [22 28 33 39 41 45 47 49 51 52 54 55];   % 55 - data is a bit noisy

nfiles=length(filenames);
col1=lines(nfiles);

case_count=0;
for jj=1:nfiles

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
  if uoo~=24 && uoo ~= 30
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

  [segments] = split_segments(hfile, uoo, deltacase);

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

  if mean_alpha>9
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
    modelalpha = model.alpha + steady_shift(index);
  else    
    modelalpha = model.alpha-1;
  end  
  modelcz = model.cz; 
 
  kall=[];
  phiall=[];
  intgall=[];
  alpha_all=[];
  theta_all=[];
  ampall=[];
  toffall=[];
  sofstall=[];
%  if (case_count>1)
%    figure(20)
%    clf
%%    close(21)
%%    close(22)
%  end

  legs_f{case_count} = [re_leg '; \alpha= ' num2str(mean_alpha), '; No:' num2str(jj)];

  col2=lines(ncases); 
 
  ksegs=0;
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

    if k<0.01
      continue
    end
   
    disp(['File: ', fname, '; No:', num2str(jj)]) 
    disp(['Reduced Frequency: ', num2str(k), ';  iseg= ', num2str(iseg)])
    disp(['Mean alpha: ', num2str(mean_alpha)])
    disp(['Amplitude: ', num2str(amp)])
    disp(['Time Period (nek): ', num2str(nek_timeperiod)])
    legs{ii} = ['Re=', num2str(Re) '; k=', num2str(k)]; 
    
    q_time = segments(iseg).qtime;
    q_alpha = segments(iseg).alpha;
  
    p_time = segments(iseg).ptime;
    p_cz = segments(iseg).Cz;
    p_cm = segments(iseg).Cm;
    
    q_cm = interp1(p_time,p_cm,q_time,'pchip');
    q_cz = interp1(p_time,p_cz,q_time,'pchip');
 
    % instantaneous rotational frequency
%    mean_alpha2 = mean_alpha*pi/180;
%    amp2 = amp*pi/180;
%    inst_alpha = mean_alpha2 + amp2*sin(omega*q_time);          % in degreees
%    inst_OMEGA = amp2*omega*cos(omega*q_time);
    
%    disp(['--------------'])
  
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
    alpha_pred = mean_alpha2 +amp2*sin(omega*q_time + theta);
    
%    figure(21)
%    plot(q_time,q_alpha*180/pi); hold on
%    plot(q_time,alpha_pred*180/pi, ' ok')
    
    % cm/cz model
    phi=-5*pi/180;
    intg_const=1;
    toff=0.;
    par0(1)=phi;
    par0(2)=intg_const;
    par0(3)=toff;
%    par0(4)=0;    % steady curve offset      
    options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-8,'TolFun', 1e-8);

%     [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model3(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

    [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model2(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

%   [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model(par,q_time,q_cz,modelcz,modelalpha,k,uoo,mean_alpha2,amp2,theta), par0,options);

%    [par,fval,exitflag,output] = fmincon(@(par) unsteady_force_model(par,q_time,q_cz,model.cz,model.alpha,k,uoo,mean_alpha2,amp2,theta),par0,[],[],[],[],[-pi/2 -100],[pi/2 100],[],options);
  
    phi=par(1);
    intg_const=par(2);
    toff=par(3);
%    steady_offset=par(4);
    
    if (isnan(fval))
      fval    
      continue
    end
    ksegs=ksegs+1;

%    [[mean_alpha2 amp2 theta phi]*180/pi intg_const]
    
    omega = 2*k*uoo/c;
    time=q_time-toff;
    
    alpha_lagg=mean_alpha2+amp2*sin(omega*time + theta + phi);
    
    inst_OMEGA=omega*amp2*cos(omega*time+theta);
    
    p_motion = intg_const*(inst_OMEGA);
    
    cz_lagg = interp1(modelalpha,modelcz,alpha_lagg*180/pi,'linear');
    
    cz_pred = p_motion + cz_lagg;
  
    kall=[kall k];
    phiall = [phiall phi];
    intgall = [intgall intg_const];
    alpha_all=[alpha_all mean_alpha2];
    theta_all=[theta_all theta];
    ampall = [ampall amp2];
    toffall = [toffall toff];    
%    sofstall = [sofstall steady_offset];  

    % Phase plot
    % interpolate onto qtime
   
    figure(10+ii)
    clf
%    subplot(ncases,1,ii)
    phase_plot = plot(q_alpha*180/pi,q_cz, ' .', 'Color', col2(ii,:)); hold on
    ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
    legend(phase_plot, legs(ii), 'Interpreter', 'tex', 'fontsize', lfs, 'Location', 'Best')
    plot(q_alpha*180/pi,cz_pred, '--', 'Color', col2(ii,:));

    if ksegs~=-1
      plot(modelalpha,modelcz,'--k', 'LineWidth', 2)
      xlim([mean_alpha2*180/pi-1.5 mean_alpha2*180/pi+1.5])
    end  

  
  end

  if found  
    figure(30)
    plot(kall,phiall*180/pi, '-.', 'Color', col1(jj,:))
    ylabel('\phi', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on    

    figure(31)
    plot(kall,intgall, '-o', 'Color', col1(jj,:))
    ylabel('Integration constant', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on
    grid on

    figure(32)
    plot(kall,toffall, '-d', 'Color', col1(jj,:))
    ylabel('Time offset', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
    hold on

%    figure(33)
%    plot(kall,sofstall, '-d', 'Color', col1(jj,:))
%    ylabel('Steady offset', 'Interpreter', 'tex', 'FontSize', fs)
%    xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
%    legend(legs_f, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
%    hold on


  end

end

%save 'all_predictions.mat'


% figure(22)
% plot(q_time,q_cz); hold on
% plot(q_time,cz_pred, ' ok')
% plot(q_time,cz_lagg, '--g')
% plot(q_time,p_motion+mean(cz_lagg), ':m')
% filename=['model_cz_time.eps'];
% filename = [re_leg '_' filename];
% %SaveFig(gcf,filename, destn, ifcols)
% 
% 
% figure(20)
% plot(alpha_pred*180/pi,cz_pred, '--m', 'LineWidth', 2)
% filename=['model_phase_potrait.eps'];
% filename = [re_leg '_' filename];
% %SaveFig(gcf,filename, destn, ifcols)




