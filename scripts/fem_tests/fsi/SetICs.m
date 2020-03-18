function [El]=SetICs(El,nelv)
% Set initial condition

[lx1 ly1]=size(El(1).xm1);
Nx=lx1-1;
Ny=ly1-1;

for elno=1:nelv
  for jj=0:Ny
    for ii=0:Nx
      xpt=El(elno).xm1(ii+1,jj+1);
      ypt=El(elno).ym1(ii+1,jj+1);
      El(elno).un(ii+1,jj+1) = usric(xpt,ypt);
      posx = jj*(Nx+1) + ii + 1;
      El(elno).unvec(posx,1) = El(elno).un(ii+1,jj+1);
    end
  end
end


