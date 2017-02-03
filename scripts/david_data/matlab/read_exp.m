% Read experimental data provided by David Eller.

%
clear
clc
close all


addpath '/scratch/negi/git_repos/matlabscripts/scripts/'

base = '/scratch/negi/exp_data/';
fol = 'delta+8/';
folder = [base fol];
lfol = length(fol)

fs = 16;

%[status,result] = system(['ls ' folder '*']);
%
%inds1 = strfind(result,fol);
%inds2 = strfind(result,'.h5');
%
%nfiles = length(inds2);
%U0=zeros(nfiles,1);
%alpha=zeros(nfiles,1)-99;
%defl=zeros(nfiles,1);
%
%disp(['N files: ', num2str(nfiles)])
%
%for i=1:nfiles
%    ind1 = inds1(i)+lfol;
%    ind2 = inds2(i)+2;
%    fname=result(ind1:ind2);
%    filenames{i}=fname;
%    inds3=strfind(fname,'_');
%    inds4=strfind(fname,'');
%   
%    if (length(inds3)>2)
%      % Freestream
%      ind1=2;
%      ind2=inds3(1)-1;
%      U0(i) = str2double(fname(ind1:ind2));
%      % apha
%      ind1=inds3(1)+2;
%      ind2=inds3(2)-1;
%      alpha(i) = str2double(fname(ind1:ind2));
%      % flap deflection      
%      ind1=inds3(2)+2;
%      ind2=inds3(3)-1;
%      defl(i) = str2double(fname(ind1:ind2));
%    else
%      % Freestream
%      ind1=2;
%      ind2=inds3(1)-1;
%      U0(i) = str2double(fname(ind1:ind2));
%      ind1=inds3(1)+1;
%      ind2=inds3(2)-1;
%      % flap deflection
%      defl(i) = str2double(fname(ind1:ind2));
%    end
%    
%end
%
%c=0.5;
%nu = 1.568E-5;


fname = ''
uoo = 24.0;
deltacase=14.0;
c=0.5;
nu = 1.568E-5;
Re=uoo*c/nu
hfile = [folder fname];

[segments] = split_segments(hfile, uoo, deltacase)

iseg = 1;
%mean_alpha = mean(segments(iseg).alpha*180/pi);
%amp = max(segments(iseg).alpha*180/pi) - mean_alpha;
%f = segments(iseg).rfreq*uoo/pi/c;
%k = segments(iseg).rfreq


nsegs = 3;

h1=figure;
h2=figure;
for iseg=1:nsegs
  figure(h1)
  plot(segments(iseg).qtime,segments(iseg).alpha*180/pi)
  hold on;
  
% plot(segments(iseg).qtime,mean_alpha + amp*cos(2*pi*f*segments(iseg).qtime), '--r')
  figure(h2)
  plot(segments(iseg).ptime,segments(iseg).Cz)
  hold on

end




