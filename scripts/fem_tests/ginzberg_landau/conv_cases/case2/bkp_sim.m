% Small test simulation.

clear
clc
close all

N=8;
periodic = 1;

[dt dt2 CONV_VAR CONV1 conv_s BC FORC abcde deriv xorg wts] = FemMat(N);

%clearvars -except N dt dt2 dc dc2 bc df_mat conv forc x abcde D

bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];

who

x=xorg;
x2 = transpose([-1:0.01:1]);

% use a normalized gaussian
%mu = 0;
%sigma = 0.2;
%u0 = normpdf(x,mu,sigma)
%u0 = u0/max(u0);
%u0_int = normpdf(x2,mu,sigma);
%u0_int = u0_int/max(u0_int);

u0 = 0.1*sin(pi*x);
u0_int = 0.1*sin(pi*x2);

%u0 = 1 - x.^2;
%u0_int = 1 - x2.^2;

un = u0;

%% testing
%xtemp = transpose(linspace(-1,1,200));
l2 = length(x2);
basis = zeros(l2,N+1);
for i = 0:N
val = zeros(l2,1);
for j = 0:N
     val = val + un(i+1)*abcde(j+1,i+1)*x2.^j;
end
     basis(:,i+1) = val;
end

vals = sum(basis,2);
h1=figure;
plot(x2,vals);
hold on
plot(x2,u0_int,'ok');
title('initial solution');

%% testing plot derivative
%xtemp = transpose(linspace(-1,1,200));
l2 = length(x2);
basis = zeros(l2,N+1);
for i = 0:N
val = zeros(l2,1);
for j = 1:N              % first index val is zero
     val = val + un(i+1)*deriv(j+1,i+1)*x2.^(j-1);
end
basis(:,i+1) = val;
end
u0_deriv = pi*0.1*cos(pi*x2);
%u0_deriv = -2*x2;

vals = sum(basis,2);
h_deriv=figure;
plot(x2,vals);
hold on
plot(x2,u0_deriv,'ok');
title('initial derivatives')

% Transformations due to periodicity
if periodic
     dt(:,1) = dt(:,1) + dt(:,end);
     dt(1,:) = dt(1,:) + dt(end,:);
     dt = dt(1:end-1,1:end-1);

     CONV1(:,1) = CONV1(:,1) + CONV1(:,end);
     CONV1(1,:) = CONV1(1,:) + CONV1(end,:);
     CONV1 = CONV1(1:end-1,1:end-1);

     dt2(:,1) = dt2(:,1) + dt2(:,end);
     dt2(1,:) = dt2(1,:) + dt2(end,:);
     dt2 = dt2(1:end-1,1:end-1);

     FORC(:,1) = FORC(:,1) + FORC(:,end);
     FORC(1,:) = FORC(1,:) + FORC(end,:);
     FORC = FORC(1:end-1,1:end-1);

     un = un(1:end-1);
     una = un;
     x_per = x(1:end-1);
end


deltat = 0.005;
istep = 0;
nsteps = 10000;

ulag1 = zeros(N+1,1);
ulag2 = zeros(N+1,1);

ualag1 = zeros(N+1,1);
ualag2 = zeros(N+1,1);


stiff = dt;
astiff = dt2;

chi = -0.0;
mass = chi*FORC - CONV1;

e1 = eig(inv(stiff)*mass);
lambda_deltat = e1*deltat
clines = load('bdfk-neutral-curve.mat');

lambdar = real(e1)*deltat;
lambdai = imag(e1)*deltat;

h2 = figure;
plot(lambdar,lambdai, '*', 'MarkerSize', 16)
hold on
plot(clines.cline3(1,2:end),clines.cline3(2,2:end), 'r')
%xlim([min(lambdar) max(lambdar)]);
%ylim([min(lambdai) max(lambdai)]);

time = 0;
xaold = x2;

%% SOLVE
h3=figure;
for i = 1:nsteps
istep = istep+1;
time = time+deltat; 

% set time integration operators
if istep == 1
     bdfk = bdf1;
     extk = ex0;

     stiff_n = stiff*bdfk(1);
     b = -(stiff*bdfk(2)*un) + deltat*extk(2)*mass*un;
     
     
     stiff_na = astiff*bdfk(1);
     ba = -(astiff*bdfk(2)*una) + deltat*extk(2)*mass*una;

elseif istep == 2
     bdfk = bdf2;
     extk = ex1;

     stiff_n = stiff*bdfk(1);
     b = -(bdfk(2)*stiff*un + bdfk(3)*stiff*ulag1) + deltat*(extk(2)*mass*un + extk(3)*mass*ulag1);
     
     stiff_na = astiff*bdfk(1);
     ba = -(bdfk(2)*astiff*una + bdfk(3)*astiff*ualag1) + deltat*(extk(2)*mass*una + extk(3)*mass*ualag1);

else
     bdfk = bdf3;
     extk = ex2;

     stiff_n = stiff*bdfk(1);
     b = -(bdfk(2)*stiff*un + bdfk(3)*stiff*ulag1 + bdfk(4)*stiff*ulag2) + deltat*(extk(2)*mass*un + extk(3)*mass*ulag1 + extk(4)*mass*ulag2);

     stiff_na = astiff*bdfk(1);
     ba = -(bdfk(2)*astiff*una + bdfk(3)*astiff*ualag1 + bdfk(4)*astiff*ualag2) + deltat*(extk(2)*mass*una + extk(3)*mass*ualag1 + extk(4)*mass*ualag2);

end

%u_extrp = 0;
%for j = 0:N
%     u_extrp = u_extrp + un(j+1)*transpose(abcde(:,j+1))*x_p_vec;
%end

un_new = gmres(stiff_n,b);
una_new = gmres(stiff_na,ba);

%una_new = stiff_na\ba;

ulag2 = ulag1;
ulag1 = un;
un = un_new;
un2 = [un; un(1)];

ualag2 = ualag1;
ualag1 = una;
una = una_new;
una2 = [una; una(1)];


%l2 = length(x2);
%basis = zeros(l2,N+1);
%for k = 0:N
%val = zeros(l2,1);
%for j = 0:N
%     val = val + un2(k+1)*abcde(j+1,k+1)*x2.^j;
%end
%basis(:,k+1) = val;
%end

%vals = sum(basis,2);
%figure(h3)
%plot(x2,vals);
%hold on
figure(h3)
plot(x,un2, 'ok', 'MarkerSize', 16);
hold on
plot(x,una2, '*r', 'MarkerSize', 16);

xanew = x2 - conv_s*time;
for j=1:length(xanew)
     if xanew(j)>1
          while(xanew(j)>1)
               xanew(j) = xanew(j)-2;
          end
     elseif xanew(j)<-1
          while (xanew(j)<-1)
               xanew(j) = xanew(j) + 2;
          end
     end
end

%u_a = 1-xanew.^2;
u_a = 0.1*sin(pi*xanew);

plot(x2,u_a,'b')
hold off
ylim([-0.15 0.15]);
pause(0.005)
%xaold = xanew;

end



