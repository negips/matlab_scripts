%%  Bunch of functions needed
%%  Local Functions

%----------------------------------------------------------------------
 
function [lsum] = GetLsum(El,varname,nelv)

  for elno=1:nelv
    evalstr = ['El(' num2str(elno) ').scrtch1 = El(' num2str(elno) ').' varname ';'];
    eval(evalstr)
  %  evalstr = ['El(' num2str(elno) ').scrtch2 = El(' num2str(elno) ').' varname ';'];
  %  eval(evalstr)
  end
  
  lsum=zeros(nelv,1);
  
  for elno=1:nelv
    lsum(elno) = sum(El(elno).scrtch1(:));
  end

return

%---------------------------------------------------------------------- 
function [l2sum] = GetL2sum(El,varname,nelv)

  for elno=1:nelv
    evalstr = ['El(' num2str(elno) ').scrtch1 = El(' num2str(elno) ').' varname ';'];
    eval(evalstr)
  %  evalstr = ['El(' num2str(elno) ').scrtch2 = El(' num2str(elno) ').' varname ';'];
  %  eval(evalstr)
  end
  
  l2sum=zeros(nelv,1);
  
  % Weighted by mass matrix gives an inner product
  for elno=1:nelv
    l2sum(elno) = sum(El(elno).scrtch1(:).*El(elno).scrtch1(:).*El(elno).mass(:));
  end

return

%---------------------------------------------------------------------- 




%%====================================================================== 
%%====================================================================== 
%% Global functions

function [glsum] = GlobSum(El,varname,nelv)

  lsum = GetLsum(El,varname,nelv);
  glsum = sum(lsum);

return

%---------------------------------------------------------------------- 
function [glnorm] = GlobL2norm(El,varname,nelv)

  l2sum = GetL2sum(El,varname,nelv);
  glnorm = sqrt(sum(l2norm));

return

%---------------------------------------------------------------------- 




%%====================================================================== 
%%====================================================================== 

%% Some Misceleneous functions

function [volume] = GetVolume(El,nelv) 

  varname='mass';
  lvol = GetLsum(El,varname,nelv);
  volume = sum(lvol);

return 
%----------------------------------------------------------------------



 
