function [lsq] = lsq_powerlaw(par,x,y)

  a=abs(par(1));
  x0=-abs(par(2));
  m=-abs(par(3));
  
  nlfit = a*(x-x0).^(m);

  lsq = mean((y-nlfit).^2);

return
