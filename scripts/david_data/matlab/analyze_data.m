% Read experimental data provided by David Eller.

%
clear
clc
close all

addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
% fname = 'u24_a3_d14_fsteps.h5';
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
    
end

c=0.5;
nu = 1.568E-5;

interesting_case_files = [];
interesting_counter = 0;

k1 = [];                % Max frequency
k2 = [];                % Second largest frequency
k1_amp = [];            % amplitude of largest frequency
k2_amp = [];            % amplitude of second largest frequency
aoa0 = [];              % mean alpha
kred_cases = [];        % Reduced frequencies
pitch_amp = [];         % pitch amplitude
re_case = [];           % Case reynolds number
re_ccode = [];          % color code for reynolds number
files_all = [];         % store all file names
hfiles_all = [];        % full path of files
seg_all = [];           % store segment index as array
U0_all = [];            % freestream
delta_all = [];         % flap deflection
allcount = 0;

red = [1 0 0];          % 573k
blue = [0 0 1];         % 765k
black = [0 0 0];        % 956k
magenta = [1 0 1];
cyan = [0 1 1];
green = [0 1 0];

for i=1:nfiles

  if alpha(i)~=-99    
    uoo = U0(i);
    deltacase=defl(1);
    Re=uoo*c/nu;
    hfile = [folder filenames{i}];

    [segments] = split_segments(hfile, uoo, deltacase);

    flds = fieldnames(segments);
    has_freq = max(strcmp(flds,'rfreq'));
    if ~has_freq
      continue
    end
%     rpm         - rpm    
%     qtime       - sampling times for alpha
%     ptime       - sampling times for pressure  (higher than qtime)
%     rfreq       - reduced frequency
%     pos         - ?
%     pressure    - pressure
%     Cz          - Normal force (ptime)
%     Cm          - Moment (ptime)?
%     delta       - flap deflection?

    rms_alpha0  = rms(segments(1).alpha - mean(segments(1).alpha));
    mean_alpha0 = mean(segments(1).alpha);
    amp0        = (rms_alpha0*sqrt(2))*180/pi;
    mean_alpha0 = mean_alpha0*180/pi;

    disp(['Filename: ', filenames{i}])
    disp(['File No: ', num2str(i)])
    disp(['Reynolds Number: ', num2str(Re)])
    disp(['Mean alpha: ', num2str(mean_alpha0)])
    disp(['********************'])

    nsegs = length(segments);
    col1 = lines(nsegs);
    
%    h1=figure;
%    h2=figure;
%    h3=figure;
%    h4=figure;
    legs = [];
    psd_all = [];
    psd_fil = [];

    tmin=1e10;
    tmax=-1;
    for iseg=1:nsegs

      allcount = allcount+1;

      rms_alpha  = rms(segments(iseg).alpha - mean(segments(iseg).alpha));
      mean_alpha = mean(segments(iseg).alpha);
      amp        = (rms_alpha*sqrt(2))*180/pi;
      mean_alpha = mean_alpha*180/pi;

      k = segments(iseg).rfreq;
      omega = 2*k*uoo/c;
      f = omega/(2*pi);
      nek_omega = 2*k;
      nek_timeperiod = 2*pi/nek_omega;

%      legs{iseg} = ['K=' num2str(k) '; iseg: ', num2str(iseg)]; 

      disp(['Reduced Frequency: ', num2str(k), ';  iseg= ', num2str(iseg)])
      disp(['Mean alpha: ', num2str(mean_alpha)])
      disp(['Amplitude: ', num2str(amp)])
      disp(['Time Period (nek): ', num2str(nek_timeperiod)])
    
      q_time = segments(iseg).qtime;
      q_alpha = segments(iseg).alpha;
%      figure(h1)
%      plot(segments(iseg).qtime,segments(iseg).alpha*180/pi, 'Color', col1(iseg,:))
%      ylabel('$\alpha$', 'Interpreter', 'Latex', 'FontSize', fs)
%      legend(legs)
%      hold on;
      
%     plot(segments(iseg).qtime,mean_alpha + amp*cos(2*pi*f*segments(iseg).qtime), '--r')
      p_time = segments(iseg).ptime;
      p_cz = segments(iseg).Cz;
      p_cm = segments(iseg).Cm;
%      figure(h2)
%      plot(p_time,p_cm, 'Color', col1(iseg,:))
%      ylabel('$C_{m}$', 'Interpreter', 'Latex', 'Fontsize', fs)
%      legend(legs)
%      hold on

      % Phase plot
      % interpolate onto qtime
      q_cz = interp1(p_time,p_cz,q_time,'pchip');
      zero_mean_q_cz = q_cz - mean(q_cz);
      norm_q_cz = zero_mean_q_cz/abs(max(zero_mean_q_cz));
      shifted_q_cz = norm_q_cz + (iseg-1)*2;

%      figure(h3)
%      plot(q_alpha*180/pi,shifted_q_cz, 'Color', col1(iseg,:))
%      ylabel('$C_{z}$', 'Interpreter', 'Latex', 'FontSize', fs)
%      xlabel('$\alpha$', 'Interpreter', 'Latex', 'FontSize', fs)
%      legend(legs)
%      hold on

      %% PSD
      tmin = min([tmin min(p_time)]);
      tmax = max([tmax max(p_time)]);

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

%      [pxx,f2] = pwelch(zero_mean_cz,window,n_overlap,nfft,p_fs);
      [pxx,f2] = pwelch(zero_mean_cm,window,n_overlap,nfft,p_fs);
      k2norm = 2*pi*f2*c/2/uoo;

      %% Find next largest frequency >1.5*fmax
      [pxx_sorted ind5] = sort(pxx,'descend');
      k2_sorted = k2norm(ind5);

      found = 0;
      ind6=1;
      k2_max = k2_sorted(1);

%     Find next largest frequency which is greater than 1.2*k2_max
%     factor is a bit arbitrary
      while ~found
        ind6=ind6+1;    
        if k2_sorted(ind6)>1.2*k2_max
           found = 1;
        end
      end

      k1 = [k1 k2_sorted(1)];                   % Max frequency
      k2 = [k2 k2_sorted(ind6)];                % Second largest frequency
      k1_amp = [k1_amp pxx_sorted(1)];          % amplitude of largest frequency
      k2_amp = [k2_amp pxx_sorted(ind6)];       % amplitude of second largest frequency
      aoa0 = [aoa0 mean_alpha];                 % mean alpha
      kred_cases = [kred_cases k];              % Reduced frequencies
      pitch_amp = [pitch_amp amp];              % pitch amplitude
      re_case = [re_case Re];                   % Reynolds number
      files_all{allcount} = filenames{i};       % filenames
      hfiles_all{allcount} = hfile;             % full path
      seg_all = [seg_all iseg];                 % segment index
      U0_all = [U0_all uoo];                    % free stream velocity
      delta_all = [delta_all deltacase];        % flap deflection 
      

      % Color code specification for Re
      if (Re<600000)
        re_ccode = [re_ccode; red];
      elseif (Re>600000 && Re < 900000)
        re_ccode = [re_ccode; blue];
      else
        re_ccode = [re_ccode; cyan];
      end        

      disp(['--------------'])
   
    end % isegs
   
  end   % if alpha~=-99

end     % fileno

ifcols = 1;
destn = 'plots';

pitch_var = max(pitch_amp) - min(pitch_amp);
dalpha_rescale = 10 + (pitch_amp - min(pitch_amp))/max(pitch_var)*10;

%%
h1=figure;
scatter(aoa0,pitch_amp,[],re_ccode, 'LineWidth', lw);
%ylim([0 2])
xlim([0 15])
ylabel('\Delta\alpha (pitch amplitude)', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
filename='aoa_pitchamp.eps'
filename = [num2str(defl(1)) '_' filename];
SaveFig(h1, filename, destn, ifcols)

%%
h2=figure;
scatter(aoa0,kred_cases,[],re_ccode, 'LineWidth', lw);
xlim([-5 15])
% errorbar(kred_cases,aoa0,pitch_amp, ' o')
ylabel('k', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
filename='k_alpha.eps'
filename = [num2str(defl(1)) '_' filename];
SaveFig(h2, filename, destn, ifcols)

%%
h3=figure;
scatter(kred_cases,k2./k1,[],re_ccode, 'LineWidth', lw);
ylabel('k2/k1', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
ylim([0 4])
grid on
filename='k_k2.eps'
filename = [num2str(defl(1)) '_' filename];
SaveFig(h3, filename, destn, ifcols)

%%
h4=figure;
scatter(kred_cases,k2_amp./k1_amp,[],re_ccode, 'LineWidth', lw);
set(gca,'Yscale', 'log')
ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
filename='k_logk2amp.eps'
filename = [num2str(defl(1)) '_' filename];
SaveFig(h4, filename, destn, ifcols)

%%
h5=figure;
scatter(kred_cases,k2_amp./k1_amp,[], re_ccode, 'LineWidth', lw);
ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
filename='k_k2amp.eps'
filename = [num2str(defl(1)) '_' filename];

SaveFig(h5, filename, destn, ifcols)

%%
h6=figure;
scatter(pitch_amp,k2_amp./k1_amp,[], re_ccode, 'LineWidth', lw);
ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
xlabel('\Delta\alpha', 'Interpreter', 'tex', 'FontSize', fs)
filename=['kratio_k2amp.eps']
filename = [num2str(defl(1)) '_' filename];
SaveFig(h6, filename, destn, ifcols)

%%
%h7=figure;
%aomega = pitch_amp.*kred_cases;
%scatter(aomega,k2_amp./k1_amp,[], re_ccode);
%ylabel('k2_{amp}/k1_{amp}', 'Interpreter', 'tex', 'FontSize', fs)
%xlabel('\alphai\omega', 'Interpreter', 'tex', 'FontSize', fs)
%filename=['aomega_k2amp.eps']
%filename = [num2str(defl(1)) '_' filename];
%SaveFig(h7, filename, destn, ifcols)

save([num2str(defl(1)) '_collated_stats.mat'])



