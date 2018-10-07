function [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_E13(e,j,il,nz,Lz,mesh2d,zf1,zf2);

   XC1 =[]; YC1 =[]; ZC1 =[]; GL1 =[]; EF1 = [];

   xt1 = zeros(4,1); yt1 = zeros(4,1); zt1 = zeros(4,1);
   xt2 = zeros(4,1); yt2 = zeros(4,1); zt2 = zeros(4,1);

   elf1 = mesh2d.cbc(1,e).connectsto;    % Element no connecting face 1
   elf3 = mesh2d.cbc(3,e).connectsto;    % Element no connecting face 3
   elf4 = mesh2d.cbc(4,e).connectsto;    % Element no connecting face 4
%  Element types
   if elf1~=0
     f1t  = mesh2d.EType{elf1};
   else
     f1t  = 'B';                   % boundary
   end
   if elf3~=0
     f3t  = mesh2d.EType{elf3};
   else
     f3t  = 'B';
   end
   if elf4~=0
     f4t  = mesh2d.EType{elf4};
   else
     f4t  = 'B';
   end  

   lz = 0;       % starting 'z'
   dz = Lz/nz;

   glno=0;

   for k=1:nz

     if mod(k,4)==1
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

     elseif mod(k,4)==2
       xt1   = mesh2d.xc(:,e);
       yt1   = mesh2d.yc(:,e);

       xt2(1) = mesh2d.xc(1,elf4);
       yt2(1) = mesh2d.yc(1,elf4);

       xt2(2) = mesh2d.xc(2,e);
       yt2(2) = mesh2d.yc(2,e);

       xt2(3) = mesh2d.xc(3,e);
       yt2(3) = mesh2d.yc(3,e);

       xt2(4) = mesh2d.xc(4,elf4);
       yt2(4) = mesh2d.yc(4,elf4);
       
       xt    = [xt1; xt2];
       yt    = [yt1; yt2];
       XC1    = [XC1 xt];
       YC1    = [YC1 yt];

       zt1 = zeros(4,1) + lz;
       lz  = lz + dz;
       zt2 = zeros(4,1) + lz;
       zt  = [zt1; zt2];
       ZC1  = [ZC1 zt];

       glno=glno+1;
       GL1 = [GL1 glno];

     elseif mod(k,4)==3

       xt1(1) = mesh2d.xc(1,elf4);
       yt1(1) = mesh2d.yc(1,elf4);

       xt1(2) = mesh2d.xc(2,e);
       yt1(2) = mesh2d.yc(2,e);

       xt1(3) = mesh2d.xc(3,e);
       yt1(3) = mesh2d.yc(3,e);

       xt1(4) = mesh2d.xc(4,elf4);
       yt1(4) = mesh2d.yc(4,elf4);

       xt2   = mesh2d.xc(:,e);
       yt2   = mesh2d.yc(:,e);

       xt    = [xt1; xt2];
       yt    = [yt1; yt2];
       XC1    = [XC1 xt];
       YC1    = [YC1 yt];

       zt1 = zeros(4,1) + lz;
       lz  = lz + dz;
       zt2 = zeros(4,1) + lz;
       zt  = [zt1; zt2];
       ZC1  = [ZC1 zt];

       glno=glno+1;
       GL1 = [GL1 glno];

     elseif mod(k,4)==0
       xt  = [mesh2d.xc(:,e); mesh2d.xc(:,e)];
       XC1  = [XC1 xt];
       yt  = [mesh2d.yc(:,e); mesh2d.yc(:,e)];
       YC1  = [YC1 yt];
       zt1 = zeros(4,1) + lz;
       lz  = lz + dz;
       zt2 = zeros(4,1) + lz;
       zt  = [zt1; zt2];
       ZC1  = [ZC1 zt];

       glno=glno+1;
       GL1 = [GL1 glno];

     end         % mod(k,4)

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

%  Build second layer


end   % function

