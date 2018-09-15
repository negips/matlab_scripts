clear
clc
close all

axfs=13;
lafs=18;

destn='plots_test/';

tr = load('tr100k.mat');

figure(7)
plot(tr.alpha, tr.trx_uv, 'b', 'LineWidth', 1.0); hold on
plot(tr.alpha, tr.trx_ww, 'r', 'LineWidth', 1.0);
smoothpvar=tr.trx_ww;
for jj=1:5
  span=10;
  smoothpvar = smooth(smoothpvar,span);
end
%ph = arrowh(tr.alpha,smoothpvar,'k',[600,90],[8 15 20 28]);
xlabel('$\alpha[^{\circ}]$')    
ylabel('$x/c$')
set(gca,'FontSize', axfs)
xlim([5.2 8.2])
legend({'$\overline{u''v''}$', '$\overline{w''w''}$'}, 'Interpreter', 'latex','FontSize', 18, 'Location', 'Best')

set(gca,'FontSize', axfs)
ylbl=get(gca,'Ylabel');
set(ylbl,'FontSize', lafs)
xlbl=get(gca,'Xlabel');
set(xlbl,'FontSize', lafs)


svfname = ['100k_transition_alpha.eps'];
SaveFig(gcf, svfname, destn, 1)

%phase lagged portrait
%phase_lag=-62.89*pi/180;
%tr.alpha2 = alpha_0 + dalpha*sin(omega*(tr.tr_time-ptch_start) + phase + phase_lag);
%
%figure(8)
%plot(tr.alpha2, tr.trx_uv, 'b', 'LineWidth', 2); hold on
% plot(tr.alpha2, tr.trx_ww, 'r', 'LineWidth', 2);
%smoothpvar=tr.trx_ww;
%for jj=1:5
%  span=10;
%  smoothpvar = smooth(smoothpvar,span);
%end
% ph = arrowh(tr.alpha2,smoothpvar,'k',[600,90],[12 17 29]);
%xlabel('$\alpha_{e}[^{\circ}]$', 'FontSize', lafs)    
%ylabel('$x/c$', 'FontSize', lafs)
%set(gca,'FontSize', axfs)
%ylbl=get(gca,'Ylabel');
%set(ylbl,'FontSize', lafs)
%xlbl=get(gca,'Xlabel');
%set(xlbl,'FontSize', lafs)
%
% legend({'$\overline{u''v''}$', '$\overline{w''w''}$'}, 'FontSize', 22, 'Location', 'Best')
%
%if ifxfoil
%  re750 = importdata('test_ed36f128+14_re7.5e5.dat');
%  
%   tr_re750_24 = interp1(re750.data(:,1),re750.data(:,7),2.4);
%   tr_re750_44 = interp1(re750.data(:,1),re750.data(:,7),4.4);
%  ind1=re750.data(:,1)>=1.5;
%  ind2=re750.data(:,1)<5;
%  ind3=find(ind1.*ind2);
%  
%  figure(8)
%  plot(re750.data(ind3,1),re750.data(ind3,7), '--k', 'LineWidth', 3); hold on
%end
%svfname = ['750k_transition_alpha_e.eps'];
%SaveFig(gcf, svfname, destn, 1)



