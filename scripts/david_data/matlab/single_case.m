% Read experimental data provided by David Eller.

%
clear
clc
close all


addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
lfol = length(fol);
destn = 'plots/';
ifcols=1;

fs = 16;
lfs = 12;

fname = 'u30_a2_d14_afsteps.h5'
%fname = 'u24_a2_d14_fsteps.h5'
uoo = 30.0;
deltacase=14.0;
c=0.5;
nu = 1.568E-5;
Re=uoo*c/nu;
hfile = [folder fname];

[segments] = split_segments(hfile, uoo, deltacase);


nsegs = length(segments);
cases=[4];
ncases=length(cases);

for ii=1:ncases

  iseg=cases(ii);
  k=segments(iseg).rfreq; 

  legs{ii} = ['K=' num2str(k) '; iseg: ', num2str(iseg)]; 

  rms_alpha  = rms(segments(iseg).alpha - mean(segments(iseg).alpha));
  mean_alpha = mean(segments(iseg).alpha);
  amp        = (rms_alpha*sqrt(2))*180/pi;
  mean_alpha = mean_alpha*180/pi;

  omega = 2*k*uoo/c;
  f = omega/(2*pi);
  nek_omega = 2*k;
  nek_timeperiod = 2*pi/nek_omega;

  disp(['Reduced Frequency: ', num2str(k), ';  iseg= ', num2str(iseg)])
  disp(['Mean alpha: ', num2str(mean_alpha)])
  disp(['Amplitude: ', num2str(amp)])
  disp(['Time Period (nek): ', num2str(nek_timeperiod)])


  figure(1)
  plot(segments(iseg).qtime,segments(iseg).alpha*180/pi)
  hold on;
  
  figure(2)
  plot(segments(iseg).ptime,segments(iseg).Cz)
  hold on

  % Phase plot
  % interpolate onto qtime
  q_time = segments(iseg).qtime;
  q_alpha = segments(iseg).alpha;
  
  p_time = segments(iseg).ptime;
  p_cz = segments(iseg).Cz;
  p_cm = segments(iseg).Cm;

  q_cz = interp1(p_time,p_cz,q_time,'pchip');
  zero_mean_q_cz = q_cz - mean(q_cz);
  norm_q_cz = zero_mean_q_cz/abs(max(zero_mean_q_cz));
  shifted_q_cz = norm_q_cz + (iseg-1)*2;

  q_cm = interp1(p_time,p_cm,q_time,'pchip');
  zero_mean_q_cm = q_cm - mean(q_cm);
  norm_q_cm = zero_mean_q_cm/abs(max(zero_mean_q_cm));
  shifted_q_cm = norm_q_cm + (iseg-1)*2;
 
  figure(3)
  if ncases~=1
    plot(q_alpha*180/pi,shifted_q_cm, 'Color', 'b')
  else
    plot(q_alpha*180/pi,q_cm, 'Color', 'b')
  end
  hold on

  figure(4)
  plot(segments(iseg).ptime,segments(iseg).Cm, 'Color', 'b')
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

  figure(5)
  pxx_norm = pxx(ind4)/max(pxx(ind4));
  pxx_shifted = pxx_norm + ii-1;
  psd_all(ii) = plot(k2(ind4),pxx_norm, 'Color', 'b'); hold on

  [val ind5] = max(pxx);    
  kmax = k2(ind5);
  ind6 = k2<1.5*kmax;
  ind7 = k2>0.5*kmax;
  ind8 = ind6.*ind7;
  pxx_fil = pxx;
  pxx_fil(find(ind8)) = 1e-15;
  pxx_fil_norm = pxx_fil/(max(pxx_fil));
  pxx_fil_shifted = pxx_fil_norm +ii-1;
  psd_fil(ii) = plot(k2(ind4),pxx_fil_norm(ind4), '--', 'Color', 'b'); hold on
  legend(psd_all, legs, 'Interpreter', 'tex', 'fontsize', fs)

  disp(['--------------'])


end

figure(1)
ylabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('time', 'Interpreter', 'tex', 'FontSize', fs)
legend(legs, 'Interpreter', 'tex', 'FontSize',lfs)
filename=['alpha_time.eps'];
filename = [num2str(deltacase) '_' filename];
SaveFig(gcf,filename, destn, ifcols)


figure(2)
ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('time', 'Interpreter', 'tex', 'FontSize', fs)
legend(legs, 'Interpreter', 'tex', 'FontSize',lfs)
filename=['cz_time.eps'];
filename = [num2str(deltacase) '_' filename];
SaveFig(gcf,filename, destn, ifcols)

figure(3)
ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
legend(legs, 'Interpreter', 'tex', 'FontSize',lfs)
filename=['cm_alpha.eps'];
filename = [num2str(deltacase) '_' filename];
SaveFig(gcf,filename, destn, ifcols)

figure(4)
ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('time', 'Interpreter', 'tex', 'FontSize', fs)
legend(legs, 'Interpreter', 'tex', 'FontSize',lfs)
filename=['cm_time.eps'];
filename = [num2str(deltacase) '_' filename];
SaveFig(gcf,filename, destn, ifcols)


figure(5)
filename=['psd.eps'];
filename = [num2str(deltacase) '_' filename];
SaveFig(gcf,filename, destn, ifcols)



