function [dnorm ddnorm] = GetDivNorm(El,nelv) 

  varname='div';
  ldiv = GetL2intg(El,varname,nelv);
  dnorm = sqrt(sum(ldiv));

  varname='divd';
  ldiv = GetDeAliasedL2intg(El,varname,nelv);
  ddnorm = sqrt(sum(ldiv));
 

return 
%----------------------------------------------------------------------

