% Load data files and plot eigenvalues
% for licenciate

clear
clc
close all

lafs = 35;
axfs = 24;

%fname='matrices_N10_NELV9_eps-2.mat';
fname{1}='matrices_N10_NELV9_eps-2_chidt1_kc1.mat';
fname{2}='matrices_N10_NELV9_eps-2_chidt10_kc1.mat';
fname{3}='matrices_N10_NELV9_eps-2_chidt50_kc1.mat';
fname{4}='matrices_N10_NELV9_eps-2_chidt100_kc1.mat';

fname{5}='matrices_N10_NELV9_eps-2_chidt1_kc3.mat';
fname{6}='matrices_N10_NELV9_eps-2_chidt10_kc3.mat';
fname{7}='matrices_N10_NELV9_eps-2_chidt50_kc3.mat';
fname{8}='matrices_N10_NELV9_eps-2_chidt100_kc3.mat';


filenames{1}='spectra_chidt_001.eps';
filenames{2}='spectra_chidt_010.eps';
filenames{3}='spectra_chidt_050.eps';
filenames{4}='spectra_chidt_100.eps';

filenames{5}='spectra_chidt_001_k8.eps';
filenames{6}='spectra_chidt_010_k8.eps';
filenames{7}='spectra_chidt_050_k8.eps';
filenames{8}='spectra_chidt_100_k8.eps';

eigfigure=figure('Units', 'normalized');
set(gcf, 'InnerPosition', [0.25 0.25 0.3 0.4])
set(gcf, 'Renderer', 'Painters')
ax1=axes;
axpos=[0.1400 0.1500 0.7650 0.7950];
set(ax1,'Position',axpos)
set(ax1,'FontSize',axfs)

cnt=0;
for ip=1:length(fname)

  cnt=cnt+1;

  load(fname{ip})
  lambda_s=diag(lambda_s);
  lambda_c=diag(lambda_c);
  
  ifbdfk=1;
  
  if ifbdfk
    lambda_c=deltat*lambda_c;
    lambda_s=deltat*lambda_s;
  end
  
  if (ifbdfk) & cnt==1 
    clines = load('bdfk-neutral-curve.mat');
    plot(clines.cline3(1,2:end),clines.cline3(2,2:end), 'k', 'LineWidth', 2); hold on
  end  
  if (ifbdfk) 
    xlabel('$\lambda_{r}\Delta t$', 'FontSize', lafs)  
    ylabel('$\lambda_{i}\Delta t$', 'FontSize', lafs)
  else
    xlabel('$\lambda_{r}$', 'FontSize', lafs)  
    ylabel('$\lambda_{i}$', 'FontSize', lafs)
  end  
  
  xlim([-1 0.1])
  ylim([-0.8 0.8])
  
  if ifbdfk
    if (cnt>1)
      delete(ipl(:))
    end    
    ipl = plot(real(lambda_s), imag(lambda_s), '.r', 'MarkerSize', 12); hold on
  else
    scatter(real(lambda_c), imag(lambda_c), '.b', 'MarkerSize', 12); hold on
    scatter(real(lambda_s), imag(lambda_s), '.r', 'MarkerSize', 12); hold on
  end
  set(gca,'LineWidth',0.75)    
  set(gca,'FontSize', axfs)
  xlbl=get(ax1,'XLabel');
  set(xlbl,'FontSize',lafs)
  ylbl=get(ax1,'YLabel');
  set(ylbl,'FontSize',lafs)
 
  %set(gca,'PlotBoxAspectRatio', [1.1 1 1])
  
  %filename=['spectra_rhs_N' num2str(Nx), '_Nxd' num2str(Nxd) '_nelv' num2str(nelv) '.eps'];
  filename=filenames{ip};
  %pause(2)
  SaveFig(eigfigure,filename,destn,1);

end
