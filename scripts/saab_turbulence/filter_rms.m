% filtering rms signals to check phase variation

clear
clc
close all

%rms = load('rms_all.mat');
load rms_all

lgfs = 8;

time = T{1};

pr = 7;
iy = 22;
s = uu{pr}(iy,:);

outpos = [0.15 0.075 0.27 0.30];
figure(1)
set(gcf,'Units', 'normalized')
set(gcf,'Position',outpos)
plot(time,s); hold on

dt = mean(diff(time));
Fs = 1/dt;
nt = length(time);

y = fft(s);

% P2 = abs(y/nt);
% P1 = P2(1:nt/2+1);
% P1(2:end-1) = 2*P1(2:end-1);

mid = nt/2+1;

nf = 790;
ind = [mid-nf:mid+nf];
yf = y;
yf(ind) = 0;

sf = ifft(yf);           % filtered signal

plot(time,sf, 'LineWidth',2)

alpha = ini_aoa + ptch_amp*sin(omega*(time-ptch_start)+ini_phase);

phase_shift = -[0 105];

outpos = [0.65 0.75 0.3 0.40];
figure(2)
set(gcf,'Units', 'normalized')
set(gcf,'Position',outpos)

for i=1:length(phase_shift)
  ps = phase_shift(i)*pi/180;
  ae = ini_aoa + ptch_amp*sin(omega*(time-ptch_start)+ini_phase+ps);
  plot(ae,sf); hold on
  lg{i}=['Phase=',num2str(phase_shift(i))];
end
legend(lg, 'Location', 'Best', 'FontSize',lgfs)
title(['Profile=',num2str(pr), ' $y_{n}=', num2str(yn{pr}(iy,1),3), '$']);

disp(['Transition lag=', num2str(-1.0976*180/pi,3)])



