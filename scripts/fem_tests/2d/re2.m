%    This has the domain and partition information
%    Right now only rectangular

% Global bounds

% Should be able to use arbitrary grid.
% But with a globally rectangular domain.

gx=zeros(4,1);
gy=zeros(4,1);

% Periodic boundaries on both sides
xperiodic = 1;
yperiodic = 1;

%%   x/y bounds for each element
%    Bottom left is element 1.
%    Right neighbor is element 2 and so on till end of row (nelx)

custom_domain=0;
if custom_domain
     sem_dd
else
     gx(1) = 0;
     gx(2) = 1;
     gx(3) = 1;
     gx(4) = 0;

     gy(1) = 0;
     gy(2) = 0;
     gy(3) = 1;
     gy(4) = 1;

     nelx=16;
     nely=16;
     nelv=nelx*nely;
     disp(['Nelx=' num2str(nelx) '; Nely=' num2str(nely) '; Nelv=' num2str(nelv)]) 

%     xvec = [0 2.5 4];
     xvec = linspace(gx(1),gx(2),nelx+1);
     yvec = linspace(gy(1),gy(4),nely+1);


     for jj=1:nely  
          for ii=1:nelx
               elno = (jj-1)*nelx + ii;
               El(elno).xc(1) = xvec((ii-1)+1);
               El(elno).xc(2) = xvec((ii-1)+2);
               El(elno).xc(3) = xvec((ii-1)+2);
               El(elno).xc(4) = xvec((ii-1)+1);

               El(elno).yc(1) = yvec((jj-1)+1);
               El(elno).yc(2) = yvec((jj-1)+1);
               El(elno).yc(3) = yvec((jj-1)+2);
               El(elno).yc(4) = yvec((jj-1)+2);

     %         Connectivity
     %          [left bottom right top]
               El(elno).conn = ElemConnectivity(elno,nelx,nely,ii,jj,xperiodic,yperiodic);
          end
     end
end

