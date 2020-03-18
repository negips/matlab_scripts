%     Testing the spectra for gibbs effect for different spatial drop offs

clear
clc
close all

n=10;
npts=2^n;
x = linspace(-10,10,npts);

X  = max(x)-min(x);     % Total period
dx = mean(diff(x));     % Sampling period
Fs = 1/dx;              % Sampling frequency                    
L  = npts;              % Length of signal


har=4;                  % harmonic
y = cos(har*pi*x);

y1 = y;
figure(1)
plot(x,y1,'LineWidth',1); hold on
ylim([-2 2])

fft1 = fft(y1);

P12 = abs(fft1/L);
P11 = P12(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);
f = Fs*(0:(L/2))/L;
omg = 2*pi*f;

figure(2)
plot(omg,P11, '.-')
xlim([0 30])
hold on


mu=1.0;
drop = exp(-(x/mu).^10);
%drop = ones(size(x));
%ind1 = abs(x)>1;
%drop(ind1) = 0;
drop = drop/max(drop);
y2 = y.*drop;
figure(1)
plot(x,y2,'LineWidth',1); hold on
plot(x,drop, '--','LineWidth',2); hold on

fft2 = fft(y2);

P22 = abs(fft2/L);
P21 = P22(1:L/2+1);
P21(2:end-1) = 2*P21(2:end-1);
%f = Fs*(0:(L/2))/L;
%omg = 2*pi*f;

figure(2)
plot(omg,P21, '.-')
%xlim([0 20])





