% Prediction of cm/cz response based on lag model

clear
clc
close all

destn = 'plots/';
ifcols= 1;

model=load('14_static_models_950k.mat');
model.alpha=model.alpha-1;

uoo=model.uoo;
deltacase=model.deltacase;
c=0.5;
nu = 1.568E-5;
Re=uoo*c/nu;
if (Re<600000)
  re_leg = '565k';
elseif (Re>600000 && Re < 900000)
  re_leg = '765k';
else
  re_leg = '950k';
end        

fs = 16;
lfs = 12;

filenames{1}='u30_a2_d14_afsteps.h5';
filenames{2}='u30_a0_d14_fsteps.h5';
filenames{3}='u30_a5_d14_afsteps.h5';
filenames{4}='u30_a7_d14_afsteps.h5';

indicies=[1];
filenames=filenames(indicies);


nfiles=length(filenames);

for jj=1:nfiles
  fname=filenames{jj};
  %fname='u24_a2_d14_fsteps.h5';
  
  hfile = ['../delta+14/' fname];
  [segments] = split_segments(hfile, uoo, deltacase);
  
  nsegs = length(segments);
  
  indicies = [1:nsegs];
  ncases = length(indicies);
  
  kall=[];
  phiall=[];
  intgall=[];

  if (jj>1)
    close(20)
    close(21)
    close(22)
  end    
 
  for ii = 1:ncases


    iseg=indicies(ii);    
    rms_alpha  = rms(segments(iseg).alpha - mean(segments(iseg).alpha));
    mean_alpha = mean(segments(iseg).alpha);
    amp        = (rms_alpha*sqrt(2))*180/pi;
    mean_alpha = mean_alpha*180/pi;
    
    k = segments(iseg).rfreq;
    omega = 2*k*uoo/c;
    f = omega/(2*pi);
    nek_omega = 2*k;
    nek_timeperiod = 2*pi/nek_omega;
    
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
    
    % Phase plot
    % interpolate onto qtime
    figure(20)
  %  subplot(ncases,1,ii)
    q_cm = interp1(p_time,p_cm,q_time,'pchip');
    q_cz = interp1(p_time,p_cz,q_time,'pchip');
   
    plot(q_alpha*180/pi,q_cz, 'Color', 'b'); hold on
    ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
    xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
    legend(legs(ii), 'Interpreter', 'tex', 'fontsize', lfs, 'Location', 'Best')
  %  hold on
  
    % instantaneous rotational frequency
    mean_alpha2 = mean_alpha*pi/180;
    amp2 = amp*pi/180;
    inst_alpha = mean_alpha2 + amp2*sin(omega*q_time);          % in degreees
    inst_OMEGA = amp2*omega*cos(omega*q_time);
  
    integrated_dp = inst_OMEGA.^2;          % times some constantant
  
    
    disp(['--------------'])
  
    theta=0;
    ini_aoa=1;
    amp2=1*pi/180;
    par0(1)=theta;
    par0(2)=ini_aoa;
    par0(3)=amp2;
    
    options=optimset('MaxFunEvals',1000,'MaxIter',10000,'TolX',1e-8);
    [par,fval,exitflag,output] = fminsearch(@(par) unsteady_alpha2(par,q_time,q_alpha,k,uoo,mean_alpha2,amp2), par0,options);
    theta=par(1);
    mean_alpha2=par(2);
    amp2=par(3);
    alpha_pred = mean_alpha2 +amp2*sin(omega*q_time + theta);
    
    figure(21)
    plot(q_time,q_alpha*180/pi); hold on
    plot(q_time,alpha_pred*180/pi, ' ok')
    
    % cm/cz model
    phi=-5*pi/180;
    intg_const=1;
    par0(1)=phi;
    par0(2)=intg_const;
    
    options=optimset('MaxFunEvals',100000,'MaxIter',100000,'TolX',1e-8,'TolFun', 1e-8);
   [par,fval,exitflag,output] = fminsearch(@(par) unsteady_force_model2(par,q_time,q_cz,model.cz,model.alpha,k,uoo,mean_alpha2,amp2,theta), par0,options);
  
%    [par,fval,exitflag,output] = fmincon(@(par) unsteady_force_model(par,q_time,q_cz,model.cz,model.alpha,k,uoo,mean_alpha2,amp2,theta),par0,[],[],[],[],[-pi/2 -100],[pi/2 100],[],options);
  
    phi=par(1);
    intg_const=par(2);
    %phi2=par(3)
    
    omega = 2*k*uoo/c;
    
    alpha_lagg=mean_alpha2+amp2*sin(omega*q_time + theta + phi);
    
    inst_OMEGA=omega*amp2*cos(omega*q_time+theta);
    
    p_motion = intg_const*(inst_OMEGA);
    
    cz_lagg = interp1(model.alpha,model.cz,alpha_lagg*180/pi,'linear');
    
    cz_pred = p_motion + cz_lagg;
  
    kall=[kall k];
    phiall = [phiall phi];
    intgall = [intgall intg_const];
  
    figure(22)
    plot(q_time,q_cz); hold on
    plot(q_time,cz_pred, ' ok')
  
  
  end
  
  figure(30)
  plot(kall,phiall*180/pi)
  hold on    


  figure(31)
  plot(kall,kall.*intgall)
  hold on

end

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




