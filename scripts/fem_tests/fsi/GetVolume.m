function [volume] = GetVolume(El,nelv) 

  varname='mass';
  lvol = GetLsum(El,varname,nelv);
  volume = sum(lvol);

return 
%----------------------------------------------------------------------

