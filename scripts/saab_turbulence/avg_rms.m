% filtering rms signals to check phase variation

clear
clc
close all

%rms = load('rms_all.mat');
load rms_all

lgfs = 8;

time = T{1};

pr = 7;

dt = mean(diff(time));
Fs = 1/dt;
nt = length(time);

mid = nt/2+1;

alpha = ini_aoa + ptch_amp*sin(omega*(time-ptch_start)+ini_phase);

% outpos = [0.15 0.075 0.27 0.30];
% figure(1)
% set(gcf,'Units', 'normalized')
% set(gcf,'Position',outpos)
% plot(time,s); hold on

tmin = min(time);
tmax = max(time);

t1 = 48.00;

tt = t1;


ind1 = time>=(t1-20*dt);
ind2 = time<=(t1+20*dt);
ind3 = find(ind1.*ind2);

nsamples = length(ind3);
disp(['No Samples: ', num2str(nsamples)])
disp(['Delta T: ', num2str(dt*nsamples)])


U1 = U{pr}(:,ind3);
U1s = mean(U1,2);

uu1  = uu{pr}(:,ind3);
uu1s = mean(uu1,2);

yn1 = yn{pr}(:,1);
gray = [0.8 0.8 0.8];

yind = 1:50;
plot(yn1(yind),uu1(yind,:), 'LineWidth', 1, 'Color', gray); hold on
plot(yn1(yind),uu1s(yind), 'LineWidth', 2, 'Color', 'k')






