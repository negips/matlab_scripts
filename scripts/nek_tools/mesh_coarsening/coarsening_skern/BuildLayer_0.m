function [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_0(e,j,il,nz,dz,Lz,mesh2d,zf1,zf2);

   XC1 =[]; YC1 =[]; ZC1 =[]; GL1 =[]; EF1 = [];

   xt1 = zeros(4,1); yt1 = zeros(4,1); zt1 = zeros(4,1);
   xt2 = zeros(4,1); yt2 = zeros(4,1); zt2 = zeros(4,1);

   lz = 0;       % starting 'z'

   glno=0;

   for k=1:nz

     xt   = [mesh2d.xc(:,e); mesh2d.xc(:,e)];
     XC1  = [XC1 xt];
     yt   = [mesh2d.yc(:,e); mesh2d.yc(:,e)];
     YC1  = [YC1 yt];
     zt1 = zeros(4,1) + lz;
     lz  = lz + dz;
     zt2 = zeros(4,1) + lz;
     zt  = [zt1; zt2];
     ZC1  = [ZC1 zt];

     glno=glno+1;
     GL1 = [GL1 glno];


     ef(1,:) = mesh2d.cbc(1,e).bc;    % BC on face 1
     ef(2,:) = mesh2d.cbc(2,e).bc;    % BC on face 2
     ef(3,:) = mesh2d.cbc(3,e).bc;    % BC on face 3
     ef(4,:) = mesh2d.cbc(4,e).bc;    % BC on face 4

     if k==1
       ef(5,:) = zf1;
       ef(6,:) = 'E  ';
     elseif k==nz;
       ef(5,:) = 'E  ';
       ef(6,:) = zf2;
     else
       ef(5,:) = 'E  ';
       ef(6,:) = 'E  ';
     end
     
     EF1{k} = ef;  

   end           % k=1:nz


end   % function

