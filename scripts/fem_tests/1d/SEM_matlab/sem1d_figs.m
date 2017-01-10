DorN = 'D'; 	% Dirichlet or Neumann
nel=60;		% number of elements
NGLL = 4:10;	% range of number of GLL nodes

if DorN =='D', e1=1; else, e1=2; end % for Neumann skip zero frequency
clf
clear dval
for ngll= NGLL,
  [eigval,err]=sem1d_homog(nel,ngll,DorN);
  %error in frequency:
  dval= abs(eigval(e1:end,1)-eigval(e1:end,2))./eigval(e1:end,2);
  %select frequencies with error > intrinsic error
  %of the eigenvalue decomposition
  isel = find(dval > 2*err./eigval(e1:end,2).^2);
  npw = 2*nel./isel; % = elements per wavelength
  loglog(npw,dval(isel)); hold on
end
hold off
title('Frequency converges as (\lambda/h)^{-2p}')
xlabel('\lambda/h = number of elements per wavelength')
ylabel('Relative error in frequency')
grid on
