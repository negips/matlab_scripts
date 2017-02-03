%     Refine stored data

clear
clc
close all

matfile = '14_collated_stats.mat';

data = load(matfile);

fs = 16;
tol_secondary= 0.05;

%Refine by k2/k1 amplitude

amp_ratio = data.k2_amp./data.k1_amp;
ind1 = find(amp_ratio>tol_secondary);        % 1%

k1r1=data.k1(ind1);
k2r1=data.k2(ind1);
k1ampr1=data.k1_amp(ind1);
k2ampr1=data.k2_amp(ind1);
aoa0r1=data.aoa0(ind1);
kredr1=data.kred_cases(ind1);
pitchr1=data.pitch_amp(ind1);
rere1=data.re_case(ind1);
filesr1=data.files_all(ind1);
segr1=data.seg_all(ind1);
reccoder1=data.re_ccode(ind1,:);
U0r1=data.U0_all(ind1);
deltar1=data.delta_all(ind1);

ifcols = 1;
destn = 'plots';

%%
h1=figure;
scatter(aoa0r1,pitchr1,[],reccoder1);
%ylim([0 2])
xlim([0 15])
ylabel('\Delta\alpha (pitch amplitude)', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
filename='aoa_pitchamp_r.eps'
filename = [num2str(deltar1(1)) '_' filename];
SaveFig(h1, filename, destn, ifcols)

%%
h2=figure;
scatter(aoa0r1,kredr1,[],reccoder1,'LineWidth', 2);
xlim([-5 15])
% errorbar(kred_cases,aoa0,pitch_amp, ' o')
ylabel('k', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
filename='k_alpha_r.eps'
filename = [num2str(deltar1(1)) '_' filename];
SaveFig(h2, filename, destn, ifcols)

%%
h3=figure;
scatter(kredr1,k2r1./k1r1,[],reccoder1);
ylabel('k2/k1', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
ylim([0 4])
grid on
filename='k_k2_r.eps'
filename = [num2str(deltar1(1)) '_' filename];
SaveFig(h3, filename, destn, ifcols)

%%
% h4=figure;
% scatter(kredr1,k2ampr1./k1ampr1,[],reccoder1);
% set(gca,'Yscale', 'log')
% ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
% xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
% filename='k_logk2amp_r.eps'
% filename = [num2str(deltar1(1)) '_' filename];
% SaveFig(h4, filename, destn, ifcols)

%%
h5=figure;
scatter(kredr1,k2ampr1./k1ampr1,[], reccoder1);
ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
filename='k_k2amp_r.eps'
filename = [num2str(deltar1(1)) '_' filename];

SaveFig(h5, filename, destn, ifcols)

%%
h6=figure;
scatter(pitchr1,k2ampr1./k1ampr1,[], reccoder1);
ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\Delta\alpha', 'Interpreter', 'tex', 'FontSize', fs)
filename=['kratio_k2amp_r.eps']
filename = [num2str(deltar1(1)) '_' filename];
SaveFig(h6, filename, destn, ifcols)

%%

save([num2str(deltar1(1)) '_collated_stats_r.mat'])

%% individual cases

amp_ratio = k2ampr1./k1ampr1;

button = 1;
counter = 0;

while button==1

  % Get case to plot
  figure(h5)
  [xp yp button] = ginput(1);
  if button~=1
    continue
  end

  % close plots
  if counter~=0
     close(20)
     close(21)
     close(22)
     close(23)
  end
  counter=counter+1;
  dist = sqrt((amp_ratio-yp).^2 + (kredr1-xp).^2);
  [val ind2] = min(dist);
  hfile = filesr1{ind2};
  iseg = segr1(ind2);
  uoo = U0r1(ind2);
  deltacase = deltar1(ind2);
  c=0.5;
  nu = 1.568E-5;
  Re=uoo*c/nu;
  col1=reccoder1(ind2,:);
% Filename      
  disp(hfile)
  [segments] = split_segments(hfile, uoo, deltacase);
  %% 

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

  legs{1} = ['\alpha=' num2str(mean_alpha) '; k=', num2str(k) '\Delta\alpha=', num2str(amp) '; Re=',num2str(Re)  ]; 

  figure(20)
  q_time = segments(iseg).qtime;
  q_alpha = segments(iseg).alpha;
  plot(segments(iseg).qtime,segments(iseg).alpha*180/pi, 'Color', col1(1,:))
  ylabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
  legend(legs, 'Interpreter', 'tex', 'fontsize', fs)
  hold on;


  figure(21)
  p_time = segments(iseg).ptime;
  p_cz = segments(iseg).Cz;
  p_cm = segments(iseg).Cm;
  plot(p_time,p_cm, 'Color', col1(1,:))
  ylabel('C_{m}', 'Interpreter', 'tex', 'Fontsize', fs)
  legend(legs, 'Interpreter', 'tex', 'fontsize', fs)
  hold on

  % Phase plot
  % interpolate onto qtime
  figure(22)
  q_cm = interp1(p_time,p_cm,q_time,'cubic');
  zero_mean_q_cm = q_cm - mean(q_cm);
  norm_q_cm = zero_mean_q_cm/abs(max(zero_mean_q_cm));
  shifted_q_cm = norm_q_cm + (iseg-1)*2;

  plot(q_alpha*180/pi,shifted_q_cm, 'Color', col1(1,:))
  ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
  xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
  legend(legs, 'Interpreter', 'tex', 'fontsize', fs)
  hold on

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

  figure(23)
  pxx_norm = pxx(ind4)/max(pxx(ind4));
  pxx_shifted = pxx_norm + iseg-1;
  psd_all = plot(k2(ind4),pxx_norm, 'Color', col1(1,:)); hold on

  [val ind5] = max(pxx);    
  kmax = k2(ind5);
  ind6 = k2<1.5*kmax;
  ind7 = k2>0.5*kmax;
  ind8 = ind6.*ind7;
  pxx_fil = pxx;
  pxx_fil(find(ind8)) = 1e-15;
  pxx_fil_norm = pxx_fil/(max(pxx_fil));
  pxx_fil_shifted = pxx_fil_norm +iseg-1;
  psd_fil = plot(k2(ind4),pxx_fil_norm(ind4), '--', 'Color', col1(1,:)); hold on
  legend(psd_all, legs, 'Interpreter', 'tex', 'fontsize', fs)

  disp(['--------------'])

end  






