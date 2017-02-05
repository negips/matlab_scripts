% Check behavior of phase lag

clear
clc
close all

fs = 16;
lfs = 12;
destn = 'plots/';
ifcols = 1;

t=linspace(0,2*pi,1000);
%phases=linspace(0,pi/2,5);
phases= [10 70 90]*pi/180;
nphase = length(phases)
col = lines(nphase);


% Phase lag with a hyperbolic tangent
% varying between [-2.5 0]
amin = -2.5;
amax = 0.;
alpha0 = (amin+amax)/2;
amp = amax - alpha0;

alpha=alpha0 + amp*sin(t);

figure(1)
alpha_all = linspace(-4,1,1000);
model = tanh(alpha_all);
plot(alpha_all,model, '--', 'LineWidth', 2); hold on
plot(alpha0,tanh(alpha0), 'ok', 'MarkerSize', 12, 'LineWidth', 2)
legend({'tanh(\alpha)'}, 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
title('Model Static Curve', 'Interpreter', 'tex', 'FontSize', lfs)
filename=['hyperbolic_tangent.eps'];
SaveFig(gcf,filename, destn, ifcols)



figure(2);
hold on;
for i=1:nphase
  phi=phases(i);
  lagalpha = alpha0 + amp*sin(t+phi);
  response = tanh(lagalpha);
  subplot(nphase,1,i)
  plot(alpha,response, 'Color', col(i,:)); hold on;
%  plot(alpha_all,model, '--')
  legs{i} = ['\phi=',num2str(phi*180/pi)];
  legend(legs(i), 'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'Best')

end
filename=['model_response.eps'];
SaveFig(gcf,filename, destn, ifcols)


%% exp_data

hfile = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/delta+14/u24_a2_d14_fsteps.h5' ;
iseg = 11;
uoo = 24;
deltacase = 14;
c=0.5;
nu = 1.568E-5;
Re=uoo*c/nu;
col1=[0 0 1];
disp(hfile)
[segments] = split_segments(hfile, uoo, deltacase);
%% 
nsegs = length(segments);

indicies = [3 7 11];
ncases = length(indicies);

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
  legs{ii} = ['k=', num2str(k)]; 
  
  q_time = segments(iseg).qtime;
  q_alpha = segments(iseg).alpha;
  
%  figure(20)
%  subplot(ncases,1,ii)
%  plot(segments(iseg).qtime,segments(iseg).alpha*180/pi, 'Color', col1(1,:))
%  ylabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
%  legend(legs(ii), 'Interpreter', 'tex', 'fontsize', fs)
%  hold on;
%  filename=['alpha_time.pdf'];
%  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)
  
  
%  figure(21)
%  subplot(ncases,1,ii)

  p_time = segments(iseg).ptime;
  p_cz = segments(iseg).Cz;
  p_cm = segments(iseg).Cm;
%  plot(p_time,p_cm, 'Color', col1(1,:))
%  ylabel('C_{m}', 'Interpreter', 'tex', 'Fontsize', fs)
%  legend(legs(ii), 'Interpreter', 'tex', 'fontsize', fs)
%  hold on
%  filename=['cm_time.pdf'];
%  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)
  
  % Phase plot
  % interpolate onto qtime
  figure(22)
  subplot(ncases,1,ii)
  q_cm = interp1(p_time,p_cm,q_time,'pchip');
%   zero_mean_q_cm = q_cm - mean(q_cm);
%   norm_q_cm = zero_mean_q_cm/abs(max(zero_mean_q_cm));
%   shifted_q_cm = norm_q_cm + (iseg-1)*2;
  
  plot(q_alpha*180/pi,q_cm, 'Color', col1(1,:))
  ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
  xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
  legend(legs(ii), 'Interpreter', 'tex', 'fontsize', lfs)
  hold on
  filename=['response_curve.eps'];
  filename = [num2str(deltacase) '_' filename];
  SaveFig(gcf,filename, destn, ifcols)
  
  %% PSD
  tmin = min(p_time);
  tmax = max(p_time);
  
  ptmax = max(segments(iseg).ptime);
  ptmin = min(segments(iseg).ptime);
  pnsamples = length(segments(iseg).ptime);
  p_sample_rate = (ptmax-ptmin)/pnsamples;
  p_fs = 1/p_sample_rate;
  disp(['Calculate frequency (p): ', num2str(p_fs)])
  p_fs = 500;
  
  qtmax = max(segments(iseg).qtime);
  qtmin = min(segments(iseg).qtime);
  qnsamples = length(segments(iseg).qtime);
  q_sample_rate = (qtmax-qtmin)/qnsamples;
  q_fs = 1/q_sample_rate;
  
  n_sections=8;        % default=8
  window = floor(pnsamples/n_sections);
  nfft = window;
  n_overlap = floor(window/2);
  zero_mean_cz = p_cz - mean(p_cz);
  zero_mean_cm = p_cm - mean(p_cm);
  
%  [pxx,f2] = pwelch(zero_mean_cz,window,n_overlap,nfft,p_fs);
  [pxx,f2] = pwelch(zero_mean_cm,window,n_overlap,nfft,p_fs);
  
  k2 = 2*pi*f2*c/2/uoo;
  
  ind4 = find(k2<2);
  
%  figure(23)
%  subplot(ncases,1,ii)
  pxx_norm = pxx(ind4)/max(pxx(ind4));
  pxx_shifted = pxx_norm + ii-1;
%  psd_all(ii) = plot(k2(ind4),pxx_norm, 'Color', col1(1,:)); hold on
  
  [val ind5] = max(pxx);    
  kmax = k2(ind5);
  ind6 = k2<1.5*kmax;
  ind7 = k2>0.5*kmax;
  ind8 = ind6.*ind7;
  pxx_fil = pxx;
  pxx_fil(find(ind8)) = 1e-15;
  pxx_fil_norm = pxx_fil/(max(pxx_fil));
  pxx_fil_shifted = pxx_fil_norm +ii-1;
%  psd_fil(ii) = plot(k2(ind4),pxx_fil_norm(ind4), '--', 'Color', col1(1,:)); hold on
%  legend(psd_all(ii), legs(ii), 'Interpreter', 'tex', 'fontsize', fs)
%  filename=['psd.pdf'];
%  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)
  
  disp(['--------------'])
end

