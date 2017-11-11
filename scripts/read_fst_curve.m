% Read turbulence intensity lineout data

clear
clc
close all

data = importdata('fst_duct_ti0001.curve');

x = data.data(:,1);
ti = data.data(:,2);

%semilogx(x, ti, 'LineWidth', 2)
plot(x, ti, 'LineWidth', 2)
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 24)
ylabel('$Ti$', 'Interpreter', 'Latex', 'FontSize', 24)
%xlim([0.1 5])
hold on

%SaveFig(gcf,'ti_semilogx.eps', 'plots/', 1)

%fop = fitoptions('Method','NonlinearLeastSquares', 'Lower', [0.001 -1 -2], 'Upper', [1 0 0], 'StartPoint', [0.05 0 -1]);
%ft=fittype('a*(x-b).^c');
%
%[f gof] = fit(x,ti, ft,'StartPoint', [0.05 0 0] )
%
%figure(2)
%plot(f,x,ti)

A=0.05;
x0=0;
m=-1;
par0(1)=A;
par0(2)=x0;
par0(3)=m;
    
options=optimset('MaxFunEvals',1000,'MaxIter',10000,'TolX',1e-8,'Tolfun',1e-8);
[par,fval,exitflag,output] = fminsearch(@(par) lsq_powerlaw(par,x,ti), par0,options)

A=abs(par(1));
x0=-abs(par(2));
m=-abs(par(3));
  
nlfit = A*(x-x0).^(m);

l1=length(nlfit);
us = [1:15:l1];         % undersample

x2=x(us);
y2=nlfit(us);

figure(1)
plot(x2,y2, 'or', 'LineWidth', 2, 'MarkerSize', 7)
lgs = ['$' num2str(round(100*A)/100) '(x+' num2str(round(100*abs(x0))/100) ')^{' num2str(round(100*m)/100) '}$'];
legend({'$Ti$', lgs}, 'FontSize', 20)


SaveFig(gcf,'ti_decay.eps', 'plots/', 1)





