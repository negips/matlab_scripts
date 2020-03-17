%% Global function
function [glsum] = GlobSum(El,varname,nelv)

  lsum = GetLsum(El,varname,nelv);
  glsum = sum(lsum);

return

%---------------------------------------------------------------------- 

