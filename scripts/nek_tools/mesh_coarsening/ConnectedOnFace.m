function newcof = ConnectedOnFace(NewCEl)

  [r nel] = size(NewCEl);
 
  newcof = []; 
  nfaces=r;
  for i=1:nel
    f=int32(zeros(nfaces,1));
    for j=1:nfaces
      if NewCEl(j,i)==0   % Boundary element
        f(j)=0;
      else
        f2=j+2;
        if f2>4
          f2=f2-4;
        end
        f(j)=f2;
      end     % if
  
    end       % j=1:nfaces
    newcof = [newcof f]; 
  
  end         % i=1:nel

end   % function  
%---------------------------------------------------------------------- 


