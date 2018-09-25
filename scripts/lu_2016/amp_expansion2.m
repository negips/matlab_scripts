% Get the growth rate from fsi_io

clear
clc
close all

lafs = 16;
lgfs = 12;

fname1 = 'fsi_io.out';
fname2 = '/scratch/negi/git_repos/fsi/stability/baseflow_solve/run_cossu/big_domain/npert_1/fsi_io.out';
fname3 = '/scratch/negi/git_repos/fsi/stability/baseflow_solve/run_cossu/big_domain/nonlinear2/fsi_io.out';

i=1; % 1
file{i}='fsi_lu_pert.out';
legen{i}='Pert';

i=i+1; % 2
file{i}='fsi_lu_base_pert.out';
legen{i}='Base+Pert';

i=i+1; % 3
file{i}='fsi_lu_base_pert_displ.out';
legen{i}='Base+Pert+Displ';

i=i+1; % 4
file{i}='fsi_lu_nl.out';
legen{i}='NL';

% Re=45
i=i+1; % 5
file{i}='fsi_lu_re45_pert.out';
legen{i}='Re45; Pert';

i=i+1; % 6
file{i}='fsi_lu_re45_base_pert.out';
legen{i}='Re45; Base+Pert';

i=i+1; % 7
file{i}='fsi_lu_re45_nl.out';
legen{i}='Re45; NL';

% After pressure bug change
i=i+1; % 8
file{i}='fsi_lu_re45_pert_2.out';
legen{i}='Re45; Pert 2';

i=i+1; % 9
file{i}='fsi_lu_re45_base_pert_2.out';
legen{i}='Re45; Base+Pert 2';

% Non-linear amplitudes

i=i+1; % 10
file{i}='fsi_lu_re45_nl_amp2.out';
legen{i}='Re45; NL; AMP2';

i=i+1; % 11
file{i}='fsi_lu_re45_nl_amp3.out';
legen{i}='Re45; NL; AMP3';

i=i+1; % 12
file{i}='fsi_lu_re45_nl_amp4.out';
legen{i}='Re45; NL; AMP4';

% Including grad(U).eta terms in forces
i=i+1; % 13
file{i}='fsi_lu_re45_all_npert2.out';
legen{i}='Re45; All; npert=2';

i=i+1; % 14
file{i}='fsi_lu_re45_all_npert1.out';
legen{i}='Re45; All; npert=1';

cols = ['b','r','k','m','c','g','y'];

% Calculate Linear growth from non-linear solutions using amplitude expansions

% dA/dt = eps1*A + eps2*A^2 + eps3*A^3

ind=[10, 11, 12];

file_in=file(ind);
legen_in=legen(ind);
nfiles = length(ind);

iskip=11;         % no of initial peaks to skip
eskip=0;          % no of end peaks to skip

%re45
tstart=1300;
tend=1400;

for i=1:nfiles

  if ind(i)==4
    eskip=8;
  end    

  fname = file_in{i};
  fsi(i) = importdata(fname);
  time{i}=fsi(i).data(:,2);
  eta{i} =fsi(i).data(:,4);
  etav{i}=fsi(i).data(:,5);

  if (tstart>0)
    ind1=time{i}>=tstart;
  else
    ind1=ones(length(time{i}),1);
  end  

  if (tend>0)
    ind2=time{i}<=tend;
  else
    ind2=ones(length(time{i}),1);
  end
  ind3=find(ind1.*ind2);
  time{i}= time{i}(ind3);
  eta{i} = eta{i}(ind3);
  etav{i} = etav{i}(ind3);
end

l=length(eta{1});
nw = 0;
ind = [];
cnd = [];

tol=1.0e-6;
tol2=1.0e-5;

for i = 2:l-1

%  if (abs(eta{1}(i))<tol && abs(eta{2}(i))<tol && abs(eta{3}(i))<tol)
%    continue
%  end

  if (abs(eta{1}(i))<tol && abs(eta{2}(i))<tol)
    continue
  end

  LHS = zeros(3,3);
  rhs = zeros(3,1);

  for j=1:3
    rhs(j,1) = (eta{j}(i+1)-eta{j}(i-1))/(time{j}(i+1)-time{j}(i-1));
    LHS(j,:) = [eta{j}(i) eta{j}(i)^2 eta{j}(i)^3];
  end 

% Derivative can be zero sometimes.   
  if (abs(rhs(1))<tol2 && abs(rhs(2))<tol2 && abs(rhs(3)<tol2))
    continue
  end

%  LHS = zeros(2,2);
%  rhs = zeros(2,1);
%
%  for j=1:2
%    rhs(j,1) = (eta{j}(i+1)-eta{j}(i-1))/(time{j}(i+1)-time{j}(i-1));
%    LHS(j,:) = [eta{j}(i) eta{j}(i)^2];
%  end

% Derivative can be zero sometimes.   
%  if (abs(rhs(1))<tol2 && abs(rhs(2))<tol2)
%    continue
%  end

  ind = [ind i];
  nw = nw + 1;

  cnd(nw) = cond(LHS);

  etav_cal(:,nw) = rhs;
  eps(:,nw) = LHS\rhs;
end

time2=time{1}(ind);
eta1 = eta{1}(ind);
etav1 = etav{1}(ind);
eps1 = eps(1,:);
eps2 = eps(2,:);

figure(1)
ax1=gca;
plot(ax1,time2,eps1, '. '); hold on
set(ax1, 'Color', 'none');
ax2=axes;
plot(ax2,time2,etav1, '.'); hold on
%plot(ax2,time2,etav_cal(1,:), 'r')
set(ax2, 'YAxisLocation', 'right');
set(ax2, 'Color', 'none');

linkaxes([ax1 ax2], 'x')
%legend(ts,legen_in,'FontSize',lgfs,'Location','Best')

figure(2)
plot(time2,cnd)

