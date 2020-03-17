function [l2intg] = GetDeAliasedL2intg(El,varname,nelv)

  EL2 = [];

  for elno=1:nelv
    evalstr = ['El2(' num2str(elno) ').scrtch1 = El(' num2str(elno) ').' varname ';'];
    eval(evalstr)
  %  evalstr = ['El(' num2str(elno) ').scrtch2 = El(' num2str(elno) ').' varname ';'];
  %  eval(evalstr)
  end
  
  l2intg=zeros(nelv,1);
  
  % Weighted by mass matrix gives an inner product
  for elno=1:nelv
    mass=diag(El(elno).massd);  
    l2intg(elno) = sum(El2(elno).scrtch1(:).*El2(elno).scrtch1(:).*mass(:));
  end

return

%---------------------------------------------------------------------- 

