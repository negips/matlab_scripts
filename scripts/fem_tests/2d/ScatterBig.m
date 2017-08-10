function [El] = ScatterBig(El,varname,bigvec,nelv)

  nbasis = length(bigvec);

  [lx1 ly1] = size(El(1).xm1);
  Nx=lx1-1;
  Ny=ly1-1;

  for elno=1:nelv
    for jj=0:Ny
      for ii=0:Nx
        gno = El(elno).lglmap(ii+1,jj+1);
        El(elno).scrtch1(ii+1,jj+1) = bigvec(gno);
      end
    end
  end

  for elno=1:nelv
    evalstr = ['El(' num2str(elno) ').' varname ' = El(' num2str(elno) ').scrtch1;'];
    eval(evalstr)
  end

  return 



 

