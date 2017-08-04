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

