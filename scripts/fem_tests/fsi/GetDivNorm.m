function [dnorm] = GetDivNorm(El,nelv) 

  varname='div';
  ldiv = GetL2intg(El,varname,nelv);
  dnorm = sqrt(sum(ldiv));

return 
%----------------------------------------------------------------------

