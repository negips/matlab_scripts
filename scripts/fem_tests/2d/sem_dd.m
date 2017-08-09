%    Custom domain decomposition for elements


gx(1) = 0;
gx(2) = 4;
gx(3) = 4;
gx(4) = 0;

gy(1) = 0;
gy(2) = 0;
gy(3) = 4;
gy(4) = 4;


nelx=2;
nely=2;
nelv=nelx*nely;

if nelx==2

%   Cartesian grid 
%     xmat(:,1) = [0; 2; 2; 0];
%     xmat(:,2) = [2; 4.0; 4.0; 2];
%     xmat(:,3) = [0; 2; 2; 0];
%     xmat(:,4) = [2; 4; 4; 2];
%
%     ymat(:,1) = [0; 0; 2; 2];
%     ymat(:,2) = [0; 0; 2; 2];
%     ymat(:,3) = [2; 2; 4; 4];
%     ymat(:,4) = [2; 2; 4; 4];


%    Approximately same size but skewed elements
     xmat(:,1) = [0; 2; 2.2; 0];
     xmat(:,2) = [2; 4.0; 4.0; 2.2];
     xmat(:,3) = [0; 2.2; 2; 0];
     xmat(:,4) = [2.2; 4; 4; 2];

     ymat(:,1) = [0; 0; 2.1; 2];
     ymat(:,2) = [0; 0; 2; 2.1];
     ymat(:,3) = [2; 2.1; 4; 4];
     ymat(:,4) = [2.1; 2; 4; 4];

%    Large Aspect ratio skewed
%     xmat(:,1) = [0; 0.3; 0.4; 0];
%     xmat(:,2) = [0.3; 4.0; 4.0; 0.4];
%     xmat(:,3) = [0; 0.4; 0.3; 0];
%     xmat(:,4) = [0.4; 4; 4; 0.3];
%
%     ymat(:,1) = [0; 0; 0.4; 0.3];
%     ymat(:,2) = [0; 0; 0.3; 0.4];
%     ymat(:,3) = [0.3; 0.4; 4; 4];
%     ymat(:,4) = [0.4; 0.3; 4; 4];

%    Large aspect ratio non-skewed
%     xmat(:,1) = [0; 0.3; 0.3; 0];
%     xmat(:,2) = [0.3; 4.0; 4.0; 0.3];
%     xmat(:,3) = [0; 0.3; 0.3; 0];
%     xmat(:,4) = [0.3; 4; 4; 0.3];
%
%     ymat(:,1) = [0; 0; 0.3; 0.3];
%     ymat(:,2) = [0; 0; 0.3; 0.3];
%     ymat(:,3) = [0.3; 0.3; 4; 4];
%     ymat(:,4) = [0.3; 0.3; 4; 4];


elseif nelx==3
     xmat(:,1) = [0; 1; 1.25; 0];
     xmat(:,2) = [1; 2; 1.8; 1.25];
     xmat(:,3) = [2; 3; 3; 1.8];
     xmat(:,4) = [0; 1.25; 1; 0];
     xmat(:,5) = [1.25; 1.8; 2; 1];
     xmat(:,6) = [1.8; 3; 3; 2];

     ymat(:,1) = [0; 0; 1.75; 2.0];
     ymat(:,2) = [0; 0; 2.1; 1.75];
     ymat(:,3) = [0; 0; 2; 2.1];
     ymat(:,4) = [2; 1.75; 4; 4];
     ymat(:,5) = [1.75; 2.1; 4; 4];
     ymat(:,6) = [2.1; 2; 4; 4];
end


[~, ny] = size(ymat);
[~, nx] = size(xmat);

if nx~=nelv || nx~=nelv
     error('nelv and custom domain inconsistent')
     break
end

for jj=1:nely  
          for ii=1:nelx
          elno = (jj-1)*nelx + ii;
          xvec = xmat(:,elno);
          yvec = ymat(:,elno);
          El(elno).xc(1) = xvec(1);
          El(elno).xc(2) = xvec(2);
          El(elno).xc(3) = xvec(3);
          El(elno).xc(4) = xvec(4);

          El(elno).yc(1) = yvec(1);
          El(elno).yc(2) = yvec(2);
          El(elno).yc(3) = yvec(3);
          El(elno).yc(4) = yvec(4);

%         Connectivity
%          [left bottom right top]
          El(elno).conn = ElemConnectivity(elno,nelx,nely,ii,jj,xperiodic,yperiodic);
     end
end
clearvars xmat ymat
