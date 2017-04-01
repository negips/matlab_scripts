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

fs = 16;
destn = 'plots/';
ifcols = 1;


[status,result] = system(['ls ' folder '*']);

inds1 = strfind(result,fol);
inds2 = strfind(result,'.h5');

nfiles = length(inds2);
U0=zeros(nfiles,1);
alpha=zeros(nfiles,1)-99;
defl=zeros(nfiles,1);

if ~isempty(strfind(fol,'+0/'))
   defl(1)=0;
elseif ~isempty(strfind(fol,'+8/'))
   defl(1)=8;   
elseif ~isempty(strfind(fol,'+11/'))
   defl(1)=11;   
elseif ~isempty(strfind(fol,'+14/'))
   defl(1)=14;   
elseif ~isempty(strfind(fol,'+24/'))
   defl(1)=24;
else
   disp('Unknow flap angle')
   break
end



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

for i=1:nfiles

  found=0;  
  if alpha(i)~=-99   && U0(i)==18
    found=1;  
    uoo = U0(i);
    deltacase=defl(1);
    Re=uoo*c/nu;
    hfile = [folder filenames{i}];

    [segments] = split_segments(hfile, uoo, deltacase);

    flds = fieldnames(segments);
    has_freq = max(strcmp(flds,'rfreq'));
    if ~has_freq
      found=0;  
      continue
    end
%     rpm         - rpm    
%     qtime       - sampling times for alpha
%     ptime       - sampling times for pressure  (higher than qtime)
%     rfreq       - reduced frequency
%     pos         - ?
%     pressure    - pressure
%     Cz          - Lift (ptime)
%     Cm          - Moment (ptime)?
%     delta       - flap deflection?
    rms_alpha0  = rms(segments(1).alpha - mean(segments(1).alpha));
    mean_alpha0 = mean(segments(1).alpha);
    amp0        = (rms_alpha0*sqrt(2))*180/pi;
    mean_alpha0 = mean_alpha0*180/pi;
    
    
    if mean_alpha0<2 || mean_alpha0>5
       found=0;
       continue
    end

    disp(['Filename: ', filenames{i}])
    disp(['File No: ', num2str(i)])
    disp(['Reynolds Number: ', num2str(Re)])
    disp(['Mean alpha: ', num2str(mean_alpha0)])
    disp(['********************'])

    nsegs = length(segments);
    col1 = lines(nsegs);
    
    h1=figure;
    h2=figure;
    h3=figure;
    h4=figure;
    legs = [];
    psd_all = [];
    psd_fil = [];

    tmin=1e10;
    tmax=-1;
    for iseg=1:nsegs

      rms_alpha  = rms(segments(iseg).alpha - mean(segments(iseg).alpha));
      mean_alpha = mean(segments(iseg).alpha);
      amp        = (rms_alpha*sqrt(2))*180/pi;
      mean_alpha = mean_alpha*180/pi;

      k = segments(iseg).rfreq;
      omega = 2*k*uoo/c;
      f = omega/(2*pi);
      nek_omega = 2*k;
      nek_timeperiod = 2*pi/nek_omega;

      legs{iseg} = ['K=' num2str(k) '; iseg: ', num2str(iseg)]; 

      disp(['Reduced Frequency: ', num2str(k), ';  iseg= ', num2str(iseg)])
      disp(['Mean alpha: ', num2str(mean_alpha)])
      disp(['Amplitude: ', num2str(amp)])
      disp(['Time Period (nek): ', num2str(nek_timeperiod)])
    
      q_time = segments(iseg).qtime;
      q_alpha = segments(iseg).alpha;
      figure(h1)
      plot(segments(iseg).qtime,segments(iseg).alpha*180/pi, 'Color', col1(iseg,:))
      ylabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
      legend(legs)
      hold on;
      
%     plot(segments(iseg).qtime,mean_alpha + amp*cos(2*pi*f*segments(iseg).qtime), '--r')
      p_time = segments(iseg).ptime;
      p_cz = segments(iseg).Cz;
      p_cm = segments(iseg).Cm;
      figure(h2)
      plot(p_time,p_cm, 'Color', col1(iseg,:))
      ylabel('C_{m}', 'Interpreter', 'tex', 'Fontsize', fs)
      legend(legs)
      hold on

      % Phase plot
      % interpolate onto qtime
      q_cz = interp1(p_time,p_cz,q_time,'pchip');
      zero_mean_q_cz = q_cz - mean(q_cz);
      norm_q_cz = zero_mean_q_cz/abs(max(zero_mean_q_cz));
      shifted_q_cz = norm_q_cz + (iseg-1)*2;

      q_cm = interp1(p_time,p_cm,q_time,'pchip');
      zero_mean_q_cm = q_cm - mean(q_cm);
      norm_q_cm = zero_mean_q_cm/abs(max(zero_mean_q_cm));
      shifted_q_cm = norm_q_cm + (iseg-1)*2;

      figure(h3)
      plot(q_alpha*180/pi,shifted_q_cm, '-', 'Color', col1(iseg,:))
      ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
      xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
      legend(legs)
      hold on

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

      k2 = 2*pi*f2*c/2/uoo;
      ind4 = find(k2<1.5);

      figure(h4)
      pxx_norm = pxx(ind4)/max(pxx(ind4));
      pxx_shifted = pxx_norm + iseg-1;
      psd_all(iseg) = plot(k2(ind4),pxx_shifted, 'Color', col1(iseg,:)); hold on
      pxx_new = pxx(ind4);
      pxx_max = max(pxx_new);

      [val ind5] = max(pxx);

      kmax = k2(ind5);
      ind6 = k2<1.5*kmax;
      ind7 = k2>0.5*kmax;
      ind8 = ind6.*ind7;
      ind9 = find(ind8);
      pxx_new = pxx;
      pxx_new(ind9) = 0;

      pxx_new_norm = pxx_new/(max(pxx_new));
      pxx_new_shifted = pxx_new_norm +iseg-1;

      psd_fil(iseg) = plot(k2(ind4),pxx_new_shifted(ind4), '--', 'Color', col1(iseg,:)); hold on
      legend(psd_all, legs, 'Location', 'SouthEast')
      title('Normalized and filtered PSD')
      xlabel('k', 'Interpreter', 'tex', 'FontSize', fs)
      ylabel('normalized psd', 'Interpreter', 'tex', 'FontSize', fs)


      disp(['--------------'])
   
    end % isegs
    figure(h1)
    plot([tmin tmax], [mean_alpha0 mean_alpha0], '--')

%     ifint = input('Save case file?: ');
%     if (ifint)
%       interesting_counter = interesting_counter+1;
%       interesting_case_files{interesting_counter} = hfile;
%     end 
   
  end   % if alpha~=-99

  if found
    figure(h2)  
    sfname='cm_time.eps';
    sfname = [num2str(deltacase) '_' sfname];
    SaveFig(h2, sfname, destn, ifcols)

    figure(h3)  
    sfname='cm_alpha.eps';
    sfname = [num2str(deltacase) '_' sfname];
    SaveFig(h3, sfname, destn, ifcols)
   
    figure(h4)  
    sfname='psd.eps';
    sfname = [num2str(deltacase) '_' sfname];
    SaveFig(h4, sfname, destn, ifcols)
  end

  validinput = 0;
  exitloop = 0;

  while (~validinput) && i~=nfiles && found==1
    inp = input('Next file(n)? or exit(e)?: ', 's')
    if strcmp(lower(inp),'n')
      validinput = 1;
      exitloop=0; 
      clc
      disp('Next file')
      close all
    elseif strcmp(lower(inp),'e')
      validinput = 1;
      exitloop=1;
    else
      validinput=0;
    end
  end

%  save('important_cases3.mat', 'interesting_case_files', '-append')
 
  if (exitloop)
    break
  end

end     % fileno


