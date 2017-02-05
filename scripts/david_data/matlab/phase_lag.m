% Check behavior of phase lag

clear
clc
close all

fs = 16;

t=linspace(0,2*pi,1000);
phases=linspace(0,pi/2,5);

nphase=length(phases);
figure(1);
hold on;
col = lines(nphase);
for i=1:nphase
  phi=phases(i);
  plot(sin(t),sin(t+phi), 'Color', col(i,:)); hold on;
  legs{i} = ['\phi=',num2str(phi*180/pi)];
end
legend(legs, 'Interpreter', 'tex', 'FontSize', fs)


% Phase lag with a hyperbolic tangent
% varying between [-2.5 0]

amin = -2.5;
amax = 0;
alpha0 = (amin+amax)/2;
amp = amax - alpha0;

alpha=alpha0 + amp*sin(t);
figure(2);
hold on;
for i=1:nphase
  phi=phases(i);
  lagalpha = alpha0 + amp*sin(t+phi);
  response = tanh(lagalpha);
  plot(alpha,response, 'Color', col(i,:)); hold on;
  legs{i} = ['\phi=',num2str(phi*180/pi)];
end
legend(legs, 'Interpreter', 'tex', 'FontSize', fs)




