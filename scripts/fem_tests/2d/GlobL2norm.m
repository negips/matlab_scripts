%% Global functions
function [glnorm] = GlobL2norm(El,varname,nelv)

  l2sum = GetL2sum(El,varname,nelv);
  glnorm = sqrt(sum(l2norm));

return

%---------------------------------------------------------------------- 

