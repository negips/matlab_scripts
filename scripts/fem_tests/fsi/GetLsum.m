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

