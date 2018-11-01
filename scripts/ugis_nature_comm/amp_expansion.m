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


ind=[10, 11, 12];

file_in=file(ind);
legen_in=legen(ind);

cols = ['b','r','k','m','c','g','y'];

iskip=11;          % no of initial peaks to skip
eskip=0;          % no of end peaks to skip
tstart=700;
tend  =920;

%re45
tstart=950;
tend=1400;

nfiles=length(ind);

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
  etav{i}= etav{i}(ind3);
 
  figure(1)
  ts(i) = plot(time{i},eta{i}, 'LineWidth', 2, 'Color', cols(i)); hold on
  xlabel('Time', 'FontSize', lafs)
  ylabel('$\eta$', 'FontSize', lafs)
  
  [pks locs] = findpeaks(eta{i});
  
  l1 = length(locs);
  locs2 = locs(iskip+1:l1-eskip);

  figure(1)  
  pks2 = pks(iskip+1:l1-eskip);
  pks_time2 = time{i}(locs2);
  ts_pks(i)=plot(pks_time2,pks2, 'o ', 'Color', cols(i), 'MarkerSize', 6); hold on
  
  tosc = diff(pks_time2);
  omg = 2*pi./tosc;
  Omega=mean(omg);
  
  growth = log(pks2(2:end)./pks2(1:end-1))./tosc;
  %growth = diff(pks2)./tosc./pks2(1:end-1);
  Growth=mean(growth);

  growth2 = log(pks2(2:end)./pks2(1:end-1));
 
  time_growth = pks_time2(2:end); 
  figure(2)
  plot(time_growth,growth, 'o-', 'Color', cols(i)); hold on
  ylabel('Growth rate')
  xlabel('Time')

  time_osc = pks_time2(2:end); 
  figure(3)
  plot(time_osc,omg, 'o-', 'Color', cols(i)); hold on
  ylabel('Angular frequency')
  xlabel('Time')

  disp(fname)
  disp(['Mean Angular frequency:', num2str(Omega,10)])
  disp(['Mean Growth rate      :', num2str(Growth,10)])

end

figure(1)
legend(ts,legen_in,'FontSize',lgfs,'Location','Best')
figure(2)
legend(legen_in,'FontSize',lgfs,'Location','Best')
figure(3)
legend(legen_in,'FontSize',lgfs,'Location','Best')

% Calculate Linear growth from non-linear solutions using amplitude expansions

% dA/dt = eps1*A + eps2*A^2 + eps3*A^3

ind=[10, 11, 12];

file_in=file(ind);
legen_in=legen(ind);

iskip=11;         % no of initial peaks to skip
eskip=0;          % no of end peaks to skip

%re45
tstart=1200;
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

end

l=length(eta{1});
nw = 0;
ind = [];

tol=1.0e-6;

for i = 2:l-1

  if (abs(eta{1}(i))<tol && abs(eta{2}(i))<tol && abs(eta{3}(i))<tol)
    continue
  end
  ind = [ind i];
  nw = nw + 1;

%  LHS = zeros(3,3);
%  rhs = zeros(3,1);
%
%  for j=1:3
%    rhs(j,1) = (eta{j}(i+1)-eta{j}(i-1))/(time{j}(i+1)-time{j}(i-1));
%    LHS(j,:) = [eta{j}(i) eta{j}(i)^2 eta{j}(i)^3];
%  end 

  LHS = zeros(2,2);
  rhs = zeros(2,1);

  for j=1:2
    rhs(j,1) = (eta{j}(i+1)-eta{j}(i-1))/(time{j}(i+1)-time{j}(i-1));
    LHS(j,:) = [eta{j}(i) eta{j}(i)^2];
  end 


  eps(:,nw) = LHS\rhs;
end



