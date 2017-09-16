% Build stability region for BDF/EXT-k

%% Third order
% beta == lambda*dt
% Polynomial == 11G^3 + (-18 -18beta)G^2 + (9 + 18beta)G + (-2 - 6beta) = 0

clear
clc
close all

npts = 1000;

betar = linspace(-1.5,0.25,npts);
betai = linspace(-1,1,npts); 

%% BDF/EXT-3
p = zeros(4,1);
for k = 1:npts
     br = betar(k);

     for j = 1:npts;
          bi = betai(j);
          z = br + 1i*bi;
          p(1) = 11;
          p(2) = -18 -18*z;
          p(3) = 9 + 18*z;
          p(4) = -2 -6*z;

          G = roots(p);

          G1(k,j) = abs(G(1));
          G2(k,j) = abs(G(2));
          G3(k,j) = abs(G(3));
          
          maxamp = max(abs(G));
          G4(k,j) = maxamp;
          if (maxamp>1)
               G5(k,j) = nan;
          else
               G5(k,j) = maxamp;
          end

     end

end
          
%h1 = figure;
%surf(betar,betai,G1, 'EdgeColor', 'none');
%view([0 90])
%colorbar
%
%h2 = figure;
%surf(betar,betai,G2, 'EdgeColor', 'none');
%view([0 90])
%colorbar
%
%h3 = figure;
%surf(betar,betai,G3, 'EdgeColor', 'none');
%view([0 90])
%colorbar

h4 = figure;
surf(betar,betai,transpose(G5), 'EdgeColor', 'none');
title('Third Order', 'FontSize', 16)
view([0 90])
colorbar

h5 = figure;
cline3 = contour(betar,betai,transpose(G4), [1 1], 'b');
title('Third Order', 'FontSize', 16)
view([0 90])
hold on
grid on

%% Second order. BDF/EXT-2
p = zeros(3,1);
for k = 1:npts
     br = betar(k);

     for j = 1:npts;
          bi = betai(j);
          z = br + 1i*bi;
          p(1) = 3;
          p(2) = -4 -4*z;
          p(3) = 1 + 2*z;

          G = roots(p);

          G1(k,j) = abs(G(1));
          G2(k,j) = abs(G(2));
          
          maxamp = max(abs(G));
          G4(k,j) = maxamp;
          if (maxamp>1)
               G5(k,j) = NaN;
          else
               G5(k,j) = maxamp;
          end

     end

end

h6 = figure;
surf(betar,betai,transpose(G5), 'EdgeColor', 'none');
title('Second Order', 'FontSize', 16)
view([0 90])
colorbar

figure(h5);
cline2 = contour(betar,betai,transpose(G4), [1 1], 'r');
title('Stability Regions', 'FontSize', 16)
view([0 90])
legend({'$3^{rd}$ order', '$2^{nd}$ Order'}, 'Interpreter', 'Latex');
%xlim([-1e-10 1e+10]);ylim([-0.1 0.1]);

clearvars -except cline2 cline3;
save('bdfk-neutral-curve.mat');

