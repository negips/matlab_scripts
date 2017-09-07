function [bigmat] = GatherBig(El,varname,nelv)

  for elno=1:nelv
    evalstr = ['El(' num2str(elno) ').scrtch1 = El(' num2str(elno) ').' varname ';'];
    eval(evalstr)
  end

  nbasis = El(1).TotDegFreedom;
  bigmat = zeros(nbasis,1);

  [lx1 ly1] = size(El(1).xm1);
  Nx=lx1-1;
  Ny=ly1-1;

  nn = 0;
  for elno=1:nelv
    for jj=0:Ny
      for ii=0:Nx
        gno = El(elno).lglmap(ii+1,jj+1);
        bigmat(gno)= bigmat(gno) + (El(elno).scrtch1(ii+1,jj+1))/El(elno).vmult(ii+1,jj+1);
      end
    end
  end

  return 



 

