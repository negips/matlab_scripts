% Parameters for time-dependent Ginzberg-Landau

clear
clc
close all


Omg = 0.1;  % Low frequency modulation
T   = 2*pi/Omg;
time = linspace(0,2*T,10000);

U0    = 2.0; 
dU    = 1.5;
U     = U0 + dU*sin(Omg*time);
gamma = 1.0;
mu_A  = (U.^2)/(4*gamma^2);

mu0   = 0.5;
dmu   = 0.25;
mu    = mu0 + dmu*cos(Omg*time); 

k_m  = mu_A/(2*gamma);

% Unstable frequencies for mean values
k = linspace(0,1,1000);
w = U0*k + 1i*(mu0 - gamma*k.^2);


disp(['U0   = ' num2str(U0,3)]);
disp(['dU   = ' num2str(dU,3)]);
disp(['Gamma= ' num2str(gamma,3)]);
disp(['mu0  = ' num2str(mu0,5)]);
disp(['dmu  = ' num2str(dmu,3)]);

disp(['Min. mu_A= ' num2str(min(mu_A),4)]);

outpos = [0.25 0.1 0.3 0.8];
figure
set(gcf,'Units','Normalized')
set(gcf,'OuterPosition',outpos)

% Peak velocity
subplot(4,1,1)
plot(time,U)
xlabel('time')
ylabel('U')

% Absolute instability k
subplot(4,1,2)
plot(time,mu_A); hold on
plot(time,mu, 'r');
legend({'$\mu_{A}$', '$\mu$'}, 'Location', 'best', 'Interpreter', 'latex')
xlabel('time')
ylabel('$\mu$')

% Maximum Growth rate
subplot(4,1,3)
plot(time,k_m);
xlabel('time')
ylabel('$k$')

% Growth rates
subplot(4,1,4)
plot(k,imag(w)); grid on
xlabel('$k$')
ylabel('$\omega_{i}$')


