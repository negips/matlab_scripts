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

ind1=re1000.data(:,1)>=alpha_min;
ind2=re1000.data(:,1)<=alpha_max;
ind3=find(ind1.*ind2);
plot(alpha_s(ind3),cz_s(ind3), 'r', 'linewidth', 4); hold on
SaveFig(gcf,static_file,'plots/',1)

omegat = linspace(0,6*pi,1000);
alpha_t = alpha0 + dalpha*sin(omegat);
cz_t = interp1(alpha_s,cz_s,alpha_t,'pchip');

figure(2)
plot(omegat/2/pi,cz_t, 'r', 'LineWidth', 2)
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









