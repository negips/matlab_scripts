% plot Cl from logfile

clear
clc
close all

fs = 16;

destn = 'plots/';   
ifcol = 1;

clfile = 'beskow_p6/cl.out';
cldata = importdata(clfile);

k=0.5;
U0=1;
semichord = 0.5;
omega = k*U0/semichord;

time = cldata.data(:,2);
alpha = 6.7 + 1.3*sin(omega*time);
area = cldata.data(:,7);
lift = cldata.data(:,4);
cl = lift./area*4;

h1=figure;
phase_cl = plot(alpha,cl, 'LineWidth', 2);
ylabel('$C_{l}$', 'Interpreter', 'Latex', 'FontSize', fs)
xlabel('$\alpha$', 'Interpreter', 'Latex', 'FontSize', fs)

svfname = 'cl_phase.eps';
SaveFig(gcf, svfname, destn, ifcol)


h2=figure;
cl_t = plot(time,cl, 'LineWidth', 2);
ylabel('$C_{l}$', 'Interpreter', 'Latex', 'FontSize', fs)
xlabel('$t$', 'Interpreter', 'Latex', 'FontSize', fs)

svfname = 'cl_time.eps';
SaveFig(gcf, svfname, destn, ifcol)

%%
cdfile = 'beskow_p6/cd.out';
cddata = importdata(cdfile);

k=0.5;
U0=1;
semichord = 0.5;
omega = k*U0/semichord;

time = cddata.data(:,2);
alpha = 6.7 + 1.3*sin(omega*time);
area = cddata.data(:,7);
drag = cddata.data(:,4);
cd = drag./area*4;

h3=figure;
phase_cd = plot(alpha,cd, 'LineWidth', 2);
ylabel('$C_{d}$', 'Interpreter', 'Latex', 'FontSize', fs)
xlabel('$\alpha$', 'Interpreter', 'Latex', 'FontSize', fs)

svfname = 'cd_phase.eps';
SaveFig(gcf, svfname, destn, ifcol)


h4=figure;
cd_t = plot(time,cd, 'LineWidth', 2);
ylabel('$C_{d}$', 'Interpreter', 'Latex', 'FontSize', fs)
xlabel('$t$', 'Interpreter', 'Latex', 'FontSize', fs)

svfname = 'cd_time.eps';
SaveFig(gcf, svfname, destn, ifcol)

