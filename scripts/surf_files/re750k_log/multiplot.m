%% Plot integral quantities

clear
clc
close all

area=0.15*1;
U0=1.0;
rho=1.0;
norm = 0.5*rho*U0^2*area;

fnames{1} = 'cl.out';
ylbls{1} = '$C_{L}$';

fnames{2} = 'cd.out';
ylbls{2} = '$C_{d}$';

fnames{3} = 'cm.out';
ylbls{3} = '$C_{m}$';

lafs = 24;        % latex fontsize

pcol = 4;         % plot variable

fileind = [1];

files = fnames(fileind);
ylbl = ylbls(fileind);
nfiles = length(files);

k=0.4;
U0=1.0;
Chord=1.0;
semichord=Chord/2;
omega=k*U0/semichord;
phase_shift=-pi/2;
ptch_start=6.0;
alpha_0=3.4;
dalpha=1.0;
Tosc=2*pi/omega;

nfigs=0;
for ii=1:nfiles

  file=files{ii};
  pvar=importdata(file);
  ind = find(pvar.data(:,1)==0);          % Remove data for istep=0.
  pvar.data(ind,:) = [];
  ind = find(pvar.data(:,2)<6.0,1,'last');  % Only plotting for after pitching start
  pvar.data(1:ind,:) = [];

  ind = find(pvar.data(:,2)>70.0);  % Only plotting for after pitching start
  pvar.data(ind,:) = [];

  nfigs=nfigs+1;
  figure(nfigs);
  ax2=axes;

  ind1=find(pvar.data(:,2)<31.5);

  plot((pvar.data(ind1,2) -ptch_start)/Tosc,pvar.data(ind1,pcol)/norm,'b', 'LineWidth', 2,'Parent', ax2); hold on
  %set(ax2, 'YLim', [1.1 1.6])

  ind2=pvar.data(:,2)>31.5;
  ind3=pvar.data(:,2)<44.85;
  ind4 = find(ind2.*ind3);
  plot((pvar.data(ind4,2)-ptch_start)/Tosc,pvar.data(ind4,pcol)/norm, 'r', 'LineWidth', 2, 'Parent',ax2)

  ind5=find(pvar.data(:,2)>=44.85);
  plot((pvar.data(ind5,2)-ptch_start)/Tosc,pvar.data(ind5,pcol)/norm, 'g', 'LineWidth', 2, 'Parent',ax2)

  set(ax2, 'Xlim', [min(pvar.data(:,2))-ptch_start max(pvar.data(:,2))-ptch_start]/Tosc)
%  set(ax2, 'Ylim', [1.1 1.42])

  ylabel(ylbl{ii}, 'FontSize', lafs, 'Parent', ax2)
  xlabel('$(t-t_{0})/T_{osc}$', 'FontSize', lafs, 'Parent', ax2)
  %set(ax2,'YColor', [0 0 1])

%ax2=axes;
  phi = omega*(pvar.data(:,2)-ptch_start);
  alpha = alpha_0 + dalpha*sin(phi+phase_shift);
  %plot(pvar.data(:,2),alpha, '--k', 'LineWidth', 1.5, 'Parent', ax2)
  %%set(ax2,'Xlim',[1 25]);
  %set(ax2,'XAxisLocation', 'top');
  %set(ax2,'YAxisLocation', 'right');
  %set(ax2,'Color', 'none')
  %set(ax2,'XColor', [1 1 1])
  %%set(ax2,'YColor', [1 0 0])
  %set(ax2, 'Xlim', [min(pvar.data(:,2)) max(pvar.data(:,2))])
  %set(ax2,'XTick', [])
  %ylabel('\alpha^{o}', 'FontSize', 20, 'Parent', ax2)
  
%  SaveFig(gcf,'cl-time-alpha750k.eps', 'plots/',1)

  nfigs=nfigs+1;    
  figure(nfigs)
  ax3=axes;
%  ind2=find(pvar.data(:,2)>0.0*pi);
  pvar2=pvar.data(ind1,pcol);
  alpha2=alpha(ind1);
  plot(alpha2,pvar2/norm, 'LineWidth',1.5, 'Parent',ax3)

  smoothpvar=pvar2/norm;
  for jj=1:5
    span=100;
    smoothpvar = smooth(smoothpvar,span);
  end
  ph = arrowh(alpha2,smoothpvar,'k',[300,60],[2 10 21 70]);
  xlabel(ax3,'$\alpha[^{\circ}]$', 'FontSize', lafs)
  ylabel(ax3,ylbl{ii}, 'FontSize', lafs)
  hold on    

%  ind3=find(pvar.data(:,2)>31.5);
  pvar3=pvar.data(ind4,pcol);
  alpha3=alpha(ind4);
  plot(alpha3,pvar3/norm, 'r', 'LineWidth', 2, 'Parent',ax3)

%  ind4=find(pvar.data(:,2)>44.85);
  pvar4=pvar.data(ind5,pcol);
  alpha4=alpha(ind5);
  plot(alpha4,pvar4/norm, 'g', 'LineWidth', 2, 'Parent',ax3)

  set(ax3, 'Ylim', [1.08 1.42])

  axpos = get(gca,'Position');    
%  set(gca,'Position', axpos + [0.02 0 -0.02 -0.02])

%  SaveFig(gcf,'cl-alpha750k.eps', 'plots/',1)

%figure(3)
%ax3=axes;
%ind=find(pvar.data(:,2)>4.0*pi);
%pvar2=pvar.data(ind,4);
%OMEGA = 1.3*omega*cos(phi);
%OMEGA2=OMEGA(ind);
%plot(OMEGA2,pvar2/norm, 'LineWidth', 1.5, 'Parent',ax3)
%xlabel('\Omega', 'FontSize', 20, 'Parent', ax3)
%ylabel('C_{L}', 'FontSize', 20, 'Parent', ax3)
%SaveFig(gcf,'pvar-omega.eps', 'plots/',1)

end
