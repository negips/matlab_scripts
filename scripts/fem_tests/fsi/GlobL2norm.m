%% Global functions
function [glnorm] = GlobL2norm(El,varname,nelv)

  l2intg = GetL2intg(El,varname,nelv);
  glnorm = sqrt(sum(l2intg));

return

%---------------------------------------------------------------------- 

