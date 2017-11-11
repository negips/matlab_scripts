% Get the slowly varying sweeps

clear
clc
close all

addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

base = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/';
fol = 'delta+14/';
folder = [base fol];
% fname = 'u24_a3_d14_fsteps.h5';
lfol = length(fol);

fs = 18;    % fontsize
lfs = 12;    % legend font size
lw = 1;     % linewidth

[status,result] = system(['ls ' folder '*alphasweep*']);

inds1 = strfind(result,fol);
inds2 = strfind(result,'.h5');

nfiles = length(inds2);
U0=zeros(nfiles,1);
defl=zeros(nfiles,1)-99;

disp(['N files: ', num2str(nfiles)])

for i=1:nfiles
    ind1 = inds1(i)+lfol;
    ind2 = inds2(i)+2;
    fname=result(ind1:ind2);
    filenames{i}=fname;

    inds3=strfind(fname,'_');
   
    if ~(isempty(inds3))
      % Freestream
      ind1=2;
      ind2=inds3(1)-1;
      U0(i) = str2double(fname(ind1:ind2));
      % flap deflection
      ind1=inds3(1)+2;
      ind2=inds3(2)-1;
      defl(i) = str2double(fname(ind1:ind2));
    end
    
end

c=0.5;
nu = 1.568E-5;

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
black = [0 0 0];        %
magenta = [1 0 1];
cyan = [0 1 1];
green = [0 1 0];        % 963k

col1 = lines(nfiles);

legs = [];

for i=2:nfiles

%  disp(['File:' filenames{i}])    
%  inp = input('Analyse file? ');
%  if ~inp
%    continue
%  end 

  ifturb = ' ';
  ind = strfind(filenames{i},'turb');
  if ~isempty(ind)
%    continue
     ifturb = 'Turb';
  end    


  if defl(i)~=-99 && U0(i) == 30
    uoo = U0(i);
    deltacase=defl(1);
    Re=uoo*c/nu;
    if (Re<600000)
      re_leg = '565k';
    elseif (Re>600000 && Re < 900000)
      re_leg = '765k';
    else
      re_leg = '950k';
    end        

    hfile = [folder filenames{i}];

    [segments] = split_segments(hfile, uoo, deltacase);
    allcount = allcount+1;
    legs{allcount} = ['Re=' re_leg '; \delta=', num2str(deltacase) '; ' ifturb];

    nsegs = length(segments);
    for iseg=1:nsegs
      qtime=segments(iseg).qtime;
      alpha=segments(iseg).alpha*180/pi;
      p_time = segments(iseg).ptime;
      p_cz = segments(iseg).Cz;
      p_cm = segments(iseg).Cm;

      figure(1)
      plot_aoa(allcount) = plot(qtime,alpha, '.', 'Color', col1(i,:)); 
      hold on

      nsmooths=5; 
      movalpha=alpha;
      for j=1:nsmooths
        order=3;
        framelen=101;  
        movalpha = sgolayfilt(movalpha,order,framelen);
      end

%      figure(2) 
%      plot(qtime,movalpha, '.', 'Color', col1(i,:));
%      hold on
   
      ptmax = max(p_time);
      ind1=find(qtime<ptmax);
      qtime2=qtime(ind1);
      alpha2=alpha(ind1);
  
      ptmin = min(p_time);
      ind2=find(qtime2>ptmin);
      qtime3=qtime2(ind2);
      alpha3=alpha2(ind2);
      q_cz = interp1(p_time,p_cz,qtime3,'pchip');
      q_cm = interp1(p_time,p_cm,qtime3,'pchip');
 
      figure(3)
      plot_cq(allcount) = plot(alpha3,q_cz, '.', 'Color', col1(i,:));
%      plot_cq(allcount) = plot(alpha3,q_cm, '.', 'Color', col1(i,:));
      ylabel('C_{z}', 'Interpreter', 'tex', 'FontSize', fs)
      xlabel('\alpha^{\circ}', 'Interpreter', 'tex', 'FontSize', fs)
      hold on
      xlim([-1 11])

      figure(4)
      plot_cq2(allcount) = plot(alpha3,q_cm, '.', 'Color', col1(i,:));
      ylabel('C_{m}', 'Interpreter', 'tex', 'FontSize', fs)
      xlabel('\alpha', 'Interpreter', 'tex', 'FontSize', fs)
      hold on
      xlim([-1 11])

    end
  end    

end    
destn = 'plots/';
ifcols = 1;

figure(1)
legend(plot_aoa,legs,'Interpreter', 'tex', 'FontSize', lfs)

figure(3)
legend(plot_cq, legs,'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'SouthEast')
grid on
filename=['cz_alpha_steady.eps'];
filename = [num2str(deltacase) '_' filename];
%SaveFig(gcf,filename, destn, ifcols)

figure(4)
legend(plot_cq2, legs,'Interpreter', 'tex', 'FontSize', lfs, 'Location', 'SouthEast')
grid on
filename=['cm_alpha_steady.eps'];
filename = [num2str(deltacase) '_' filename];
%SaveFig(gcf,filename, destn, ifcols)


ifunsteady = 0;

if (ifunsteady)

  hfile = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/delta+14/u24_a2_d14_fsteps.h5';
  uoo = 24;
  deltacase=14;
  iseg=11;
        
  [segments] = split_segments(hfile, uoo, deltacase);

  q_time = segments(iseg).qtime;
  q_alpha = segments(iseg).alpha;

  p_time = segments(iseg).ptime;
  p_cm = segments(iseg).Cm;
  p_cz = segments(iseg).Cz;

  figure(3)
  q_cm = interp1(p_time,p_cm,q_time,'pchip');
  q_cz = interp1(p_time,p_cz,q_time,'pchip');
  %  zero_mean_q_cm = q_cm - mean(q_cm);
%  norm_q_cm = zero_mean_q_cm/abs(max(zero_mean_q_cm));
%  shifted_q_cm = norm_q_cm + (iseg-1)*2;
  plot(q_alpha*180/pi,q_cz, 'Color', 'k')
  filename=['cz_alpha_steady_unsteady.eps'];
  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)

  figure(4)
  %  zero_mean_q_cm = q_cm - mean(q_cm);
%  norm_q_cm = zero_mean_q_cm/abs(max(zero_mean_q_cm));
%  shifted_q_cm = norm_q_cm + (iseg-1)*2;
  plot(q_alpha*180/pi,q_cm, 'Color', 'k')
  filename=['cm_alpha_steady_unsteady.eps'];
  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)


end

%% David Eller's suggested case
ifdavid = 0;

if (ifdavid)

  hfile = '/scratch/negi/git_repos/matlabscripts/scripts/david_data/delta+14/u30_a2_d14_fsteps.h5';
  uoo = 30;
  deltacase=14;
  iseg=4;
        
  [segments] = split_segments(hfile, uoo, deltacase);

  q_time = segments(iseg).qtime;
  q_alpha = segments(iseg).alpha;

  p_time = segments(iseg).ptime;
  p_cm = segments(iseg).Cm;
  p_cz = segments(iseg).Cz;

  figure(3)
  q_cz = interp1(p_time,p_cz,q_time,'pchip');
  q_cm = interp1(p_time,p_cm,q_time,'pchip');
  plot(q_alpha*180/pi,q_cz, 'Color', 'k')
  filename=['cz_alpha_steady_david.eps'];
  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)

  figure(4)
  q_cz = interp1(p_time,p_cz,q_time,'pchip');
  q_cm = interp1(p_time,p_cm,q_time,'pchip');
  plot(q_alpha*180/pi,q_cm, 'Color', 'k')
  filename=['cm_alpha_steady_david.eps'];
  filename = [num2str(deltacase) '_' filename];
%  SaveFig(gcf,filename, destn, ifcols)
end



%%
% Xfoil data

ifxfoil = 0;

if ifxfoil
  xfile = 'polar_re1e6_ed36f128+14.dat';
  
  xfoil = importdata(xfile);
  figure(3)
  plot(xfoil.data(:,1), xfoil.data(:,2), '--k', 'LineWidth', 3)

  figure(4)
  plot(xfoil.data(:,1), xfoil.data(:,5), '--k', 'LineWidth', 3)

end



