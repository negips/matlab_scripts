%% Damping function used for mesh rotation

clear
clc
close all


lafs=18;

h1 = linspace(0,2,1000);
h2 = h1;
y2=-2;
y1=2.0;
x1=0.35;
x2=1.75;
ind=find(h1<x1);
h2(ind)=x1;
ind=find(h1>x2);
h2(ind)=x2;

rr = (y2-y1)/(x2-x1).*h2 + (y1*x2-y2*x1)/(x2-x1);

f = (erf(rr)-erf(y2))/(erf(y1) - erf(y2));

plot(h1,f, 'LineWidth', 2)
xlabel('$|r|$', 'FontSize', lafs, 'Interpreter', 'Latex')
ylabel('$f$', 'FontSize', lafs, 'Interpreter', 'Latex')
SaveFig(gcf,'damping_func.eps', 'plots/',1)



