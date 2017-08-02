function [evec lambda] = SystemEig(sysmat,ifplot,plothandle,ifsparsity,sparsityhandle,ifstability,col)


  destn = './plots';

  if (ifsparsity) 
    figure(sparsityhandle)
    spytol(sysmat,6,col)

    % filename = 'spy_dealiased_N4_nelv4_uniform';
    % SaveFig(h, filename, destn, 1)
  end

%% Eigenvalues:
  deltat = 1.0;
  [evec lambda] = eig(sysmat);

  lambdar = real(lambda)*deltat;
  lambdai = imag(lambda)*deltat;

  if ifplot
    figure(plothandle);
    plot(lambdar,lambdai,'.', 'Color',col, 'MarkerSize', 8)
    hold on
    if (ifstability)  
      clines = load('bdfk-neutral-curve.mat');
      plot(clines.cline3(1,2:end),clines.cline3(2,2:end), 'r')
    end  

    %xlim([min(lambdar) max(lambdar)]);
    %ylim([min(lambdai) max(lambdai)]);
  end


