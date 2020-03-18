% Build a static response in time.

clear
clc
close all

lafs=24;
axfs=13;
iflinear=0;

if iflinear
  static_file='linear_cz_re1e6.eps';
  dynamic_file='dynamic_linear_cz_re1e6.eps';
else
  static_file='nonlinear_cz_re1e6.eps';
  dynamic_file='dynamic_nonlinear_cz_re1e6.eps';
end

re100 = importdata('test_ed36f128+14_re1e5.dat');
re750 = importdata('test_ed36f128+14_re7.5e5.dat');
re1000 = importdata('test_ed36f128+14_re1e6.dat');

alpha_s = re1000.data(:,1);   % static
cz_s    = re1000.data(:,2);   % static

figure(1)
set(gcf,'Units','normalized')
pos=[0.05 0.20 0.30 0.39];
set(gcf,'Position',pos)
%plot(re100.data(:,1),re100.data(:,2), 'b', 'linewidth', 1.5); hold on
%plot(re750.data(:,1),re750.data(:,2), 'r', 'linewidth', 1.5);
plot(re1000.data(:,1),re1000.data(:,2), 'k', 'linewidth', 2); hold on
xlim([0 7])
%grid on
ylabel('$C_{z}$', 'FontSize', lafs)
xlabel('$\alpha[^{\circ}]$', 'FontSize', lafs)
set(gca,'FontSize', axfs)
SaveFig(gcf,'static_cz_re1e6.eps','plots/',1)

%x1=0.5;
%x2=2.5;
%
%ylims = get(gca,'YLim');
%X = [x1 x1 x2 x2];
%Y = [ylims fliplr(ylims)];
%f1 = fill(X,Y,'b');
%ylim(ylims)
%alpha(f1,0.1)


if iflinear
  alpha0 = 1.5;
  dalpha = 1.0;
else
  alpha0 = 2.7;
  dalpha = 1.0;
end

alpha_min = alpha0 - dalpha;
alpha_max = alpha0 + dalpha;

ind1=re1000.data(:,1)>=alpha_min-0.001;
ind2=re1000.data(:,1)<=alpha_max;
ind3=find(ind1.*ind2);
plot(alpha_s(ind3),cz_s(ind3), 'r', 'linewidth', 4); hold on
SaveFig(gcf,static_file,'plots/',1)


nt=1000;
omegat = linspace(0,6*pi,nt);
alpha_t = alpha0 + dalpha*sin(omegat);
cz_t = interp1(alpha_s,cz_s,alpha_t,'pchip');

figure(2)
set(gcf,'Units','normalized')
pos=[0.40 0.20 0.30 0.39];
set(gcf,'Position',pos)
plot(omegat/2/pi,cz_t, 'r', 'LineWidth', 2); hold on
if iflinear
  ylim([1.0 1.3])
else
  ylim([1.15 1.30])
end
ylb = ylabel('$C_{z}$', 'FontSize', lafs);
xlb = xlabel('$t/T_{osc}$', 'FontSize', lafs);
set(gca,'FontSize', axfs)
set(xlb,'Units', 'Normalized')
xpos = get(xlb,'Position');
set(xlb,'Position',xpos + [0 0.02 0])
SaveFig(gcf,dynamic_file,'plots/',1)


% Time dependent
tstmp = omegat/2/pi;

figure(3)
set(gcf,'Units','normalized')
pos=[0.20 0.20 0.90 0.60];
set(gcf,'Position',pos)

sp1 = subplot(1,3,1);
set(sp1,'Position', [0.05 0.11 0.20 0.8150])
set(sp1,'Box','on');
hold on

amp = 0.015;
theta = -50/180*pi;
ell = amp*sin(omegat+theta);
sp1_xlims = [min(alpha_t) max(alpha_t)];
sp1_ylims = 1.1*[-0.5 0.5];

sp2 = subplot(1,3,2);
set(sp2,'Position', [0.35 0.11 0.20 0.8150])
set(sp2,'Box','on');
plot(re1000.data(:,1),re1000.data(:,2), 'k', 'linewidth',1); hold on
plot(alpha_s(ind3),cz_s(ind3), 'r', 'linewidth', 1.5); hold on
xlim([1 6])
set(sp2,'XTick',[])
set(sp2,'YTick',[])

sp3 = subplot(1,3,3);
set(sp3,'Position', [0.65 0.11 0.30 0.8150])
set(sp3,'Box','on');
hold on

czz = ell + cz_t;

dim=[0.59 0.48 0.1 0.1];
str = '=';
ann1 = annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',24,'LineStyle','none');

dim=[0.29 0.48 0.1 0.1];
str = '+';
ann1 = annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',24,'LineStyle','none');


for i=1:nt
  if (i>1)

    delete(dpl1)
    delete(dpp1)
        
    delete(dpl3)
    delete(dpp3)
    
    delete(dpp2)
  end

  axes(sp1)
  dpl1 = plot(alpha_t(1:i),ell(1:i),'-r','LineWidth',1.5); 
  dpp1 = plot(alpha_t(i),ell(i),'or','LineWidth',1.5);

  xlim(sp1,sp1_xlims)
  ylim(sp1,sp1_ylims)
  set(sp1,'XTick',[])
  set(sp1,'YTick',[])

  dpl3 = plot(sp3,tstmp(1:i),czz(1:i),'-r','LineWidth',1.5);
  dpp3 = plot(sp3,tstmp(i),czz(i),'or','LineWidth',1.5);
  if iflinear
    ylim(sp3,[1.0 1.3])
  else
    ylim(sp3,[1.15-amp 1.30+amp])
  end
  xlim(sp3,[0 3])
  set(sp3,'XTick',[])
  set(sp3,'YTick',[])

  dpp2 = plot(sp2,alpha_t(i),cz_t(i),'or','LineWidth',2);

  pause(0.001)
end








