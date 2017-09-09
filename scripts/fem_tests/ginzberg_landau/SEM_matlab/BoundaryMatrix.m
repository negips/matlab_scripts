function [B,iB] = BoundaryMatrix(wgll,NELXY,iglob,jac1D,side)

NGLL = length(wgll);
NELX = NELXY(1);
NELY = NELXY(2);

switch side

 % left
  case 'L'
    ng = NELY*(NGLL-1)+1;
    iB = zeros(ng,1);
    B  = zeros(ng,1);
    for ey=1:NELY,
      ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
      e=(ey-1)*NELX+1;
      iB(ip) = iglob(1,1:NGLL,e);
      B(ip) = B(ip) + jac1D*wgll;
    end

  case 'R'
    ng = NELY*(NGLL-1)+1;
    iB = zeros(ng,1);
    B = zeros(ng,1);
    for ey=1:NELY,
      ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
      e=(ey-1)*NELX+NELX;
      iB(ip) = iglob(NGLL,1:NGLL,e);
      B(ip) = B(ip) + jac1D*wgll;
    end

  case 'T'
    ng = NELX*(NGLL-1)+1;
    iB = zeros(ng,1);
    B = zeros(ng,1);
    for ex=1:NELX,
      ip = (NGLL-1)*(ex-1)+[1:NGLL] ;
      e=(NELY-1)*NELX+ex;
      iB(ip) = iglob(1:NGLL,NGLL,e);
      B(ip) = B(ip) + jac1D*wgll;
    end

  case 'B'
    ng = NELX*(NGLL-1)+1;
    iB = zeros(ng, 1); 
    B = zeros(ng, 1); 
    for ex=1:NELX,
      ip = (NGLL-1)*(ex-1)+[1:NGLL];
      e = ex;
      iB(ip) = iglob(1:NGLL,1,e);
      B(ip) = B(ip) + jac1D*wgll;
    end

end
