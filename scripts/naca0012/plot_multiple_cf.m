%     Load a sved surface and plot individual cf

clc
close all

%load('re100k_surface.mat');


%figure(1)
%plot(surf_t11(:,1));

tnorm = (surf_t11(:,1)-ptch_start)/Tosc;

ti = [3.27 3.32 3.4];

cols = lines(length(ti));

figure(2)
for it=1:length(ti)
  [val ind] = min(abs(tnorm - ti(it)));
  plot(surf_x11(ind,:),surf_v11(ind,:), '.', 'Color', cols(it,:)); hold on
end

