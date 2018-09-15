%% Plot some individual phases

clear
clc

load('re750k_surface.mat')

tin = ([1:4]+0.6)*Tosc;
osc_start = 6.0;

close all

cols = lines(length(tin));

for i=1:length(tin)

  [t5diff ind5] = min(abs(surf_t5(:,1) - osc_start - tin(i)));
  [t8diff ind8] = min(abs(surf_t8(:,1) - osc_start - tin(i)));

  if t5diff<t8diff
    phase_t = surf_t5(ind5,:);
    phase_x = surf_x5(ind5,:);
    phase_v = surf_v5(ind5,:);
    phase_p = surf_p5(ind5,:);
  else
    phase_t = surf_t8(ind8,:);
    phase_x = surf_x8(ind8,:);
    phase_v = surf_v8(ind8,:);
    phase_p = surf_p8(ind8,:);
  end

  legs{i} = ['T=' num2str((phase_t(1)-osc_start)/Tosc)];   
      
  figure(1)
  plot(phase_x,phase_p, '-',  'Color', cols(i,:), 'MarkerSize', 12)
  hold on
end
legend(legs, 'FontSize', 16, 'Location', 'Best')
xlim([0 1]) 
