%     Analyze the fsi_terms.out file

clear
clc
close all

fsi = importdata('fsi_terms.out');
eta = importdata('fsi_io.out');

dt = 0.8e-4;
ist = fsi.data(:,1);
fs = fsi.data(:,2);
fg = fsi.data(:,3);
alpha = fsi.data(:,5);

ft = fs + alpha.*fg;



ind = 1000001:2000000;

figure(1)
plot(ist(ind),ft(ind))


f = 1/dt;             % Sampling frequency                    
L = length(ind);      % Length of signal
T = dt;               % Sampling period       
t = (0:L-1)*T;        % Time vector

Y = fft(ft(ind));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

freqs = f*(0:(L/2))/L;

figure(2)
loglog(freqs,P1) 
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
ylabel('$|P1(f)|$')

Kred = 0.183;
Omega = Kred*1.0/0.5;

% positive frequencies
N = length(ind);
Yp = Y(2:N/2+1);
Yn = Y(N/2+2:end);


fcut = 3.0*Omega/(2*pi);
index = find(freqs>fcut,1);

r = index;

% r = 20; % range of frequencies we want to preserve

mask = zeros(length(Y),1);
mask(1:r+1) = 1;
mask(end:-1:end-r+1)=1;

Yfil = Y.*mask;

y2 = ifft(Yfil);
figure(1)
hold on
plot(ist(ind),y2, 'LineWidth', 2)


alpha = 180/pi*eta.data(ind,4);
time = eta.data(ind,2);
[pk loc] = findpeaks(alpha);
tpks = time(loc);
Tosc = mean(diff(tpks));
Omega=2*pi/Tosc;

figure(3)
plot(alpha,y2)
phi=102*pi/180;
alpha2 = 5.5*sin(Omega*time + phi);

phi2 = 10*pi/180;
amp = -0.02*cos(Omega*time + phi + phi2);

figure(4)
plot(alpha2,y2); hold on
plot(alpha2,y2-amp, 'r')




