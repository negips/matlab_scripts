function [XC1,YC1,ZC1,GL1,j]=BuildLayer_ER(XC1,YC1,ZC1,GL1,e,j,il,nz,Lz,mesh2d)

%  Only one layer needs to be built
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
   nz=nz/2;

   for k=1:nz

     lz1=lz+dz;
     lz2=lz+2*dz;

     if mod(k,2)==1
       xt  = [mesh2d.xc(:,e); mesh2d.xc(:,e)];
       XC1 = [XC1 xt];
       yt  = [mesh2d.yc(:,e); mesh2d.yc(:,e)];
       YC1 = [YC1 yt];

       zt1 = zeros(4,1) + lz;

       zt2(1) = lz2;
       zt2(2) = lz2;
       zt2(3) = lz1;
       zt2(4) = lz1;

       zt  = [zt1; zt2];
       ZC1 = [ZC1 zt];

       glno=glno+1;
       GL1 = [GL1 glno];

       lz = lz2;

     else
       xt  = [mesh2d.xc(:,e); mesh2d.xc(:,e)];
       XC1 = [XC1 xt];
       yt  = [mesh2d.yc(:,e); mesh2d.yc(:,e)];
       YC1 = [YC1 yt];

       zt1(1) = lz;
       zt1(2) = lz;
       zt1(3) = lz1;
       zt1(4) = lz1;

       zt1 = zeros(4,1) + lz2;

       zt  = [zt1; zt2];
       ZC1 = [ZC1 zt];

       glno=glno+1;
       GL1 = [GL1 glno];

       lz = lz2;

     end         % if mod(k,2)

   end           % k=1:nz

end   % function

