function [XC1,XC2,YC1,YC2,ZC1,ZC2,GL1,GL2]=BuildLayer_S(XC1,XC2,YC1,YC2,ZC1,ZC2,GL1,GL2,e,j,il,nz,Lz,mesh2d);

%  Build the first layer
   xt1 = zeros(4,1); yt1 = zeros(4,1); zt1 = zeros(4,1);
   xt2 = zeros(4,1); yt2 = zeros(4,1); zt2 = zeros(4,1);
   xt3 = zeros(4,1); yt3 = zeros(4,1); zt3 = zeros(4,1);
   xt4 = zeros(4,1); yt4 = zeros(4,1); zt4 = zeros(4,1);

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

       if strcmpi(f1t,'e4')
         xt2(2) = mesh2d.xc(2,elf1);
         yt2(2) = mesh2d.yc(2,elf1);
       else
         xt2(2) = mesh2d.xc(2,e);
         yt2(2) = mesh2d.yc(2,e);
       end

       if strcmpi(f3t,'e4')
         xt2(3) = mesh2d.xc(3,elf3);
         yt2(3) = mesh2d.yc(3,elf3);
       else
         xt2(3) = mesh2d.xc(3,e);
         yt2(3) = mesh2d.yc(3,e);
       end
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
       if strcmpi(f1t,'e4')
         xt1(2) = mesh2d.xc(2,elf1);
         yt1(2) = mesh2d.yc(2,elf1);
       else
         xt1(2) = mesh2d.xc(2,e);
         yt1(2) = mesh2d.yc(2,e);
       end

       if strcmpi(f3t,'e4')
         xt1(3) = mesh2d.xc(3,elf3);
         yt1(3) = mesh2d.yc(3,elf3);
       else
         xt1(3) = mesh2d.xc(3,e);
         yt1(3) = mesh2d.yc(3,e);
       end
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

   end           % k=1:nz


%  Build second layer

   nz=nz/2;
   dz2=dz*2;
   lz=0;
   
   e=elf4;     

   for k=1:nz

     if mod(k,2)==1
       xt1  = mesh2d.xc(:,e);
       yt1  = mesh2d.yc(:,e);

       xt2  = xt1;
       yt2  = yt1;

       xt = [xt1; xt2];
       yt = [yt1; yt2];

       XC2  = [XC2 xt];
       YC2  = [YC2 yt];

       zt1    = zeros(4,1) + lz;
       lz2    = lz + dz2;
       lz1    = lz + dz;
       zt2(1) = lz2;
       zt2(2) = lz1;
       zt2(3) = lz1;
       zt2(4) = lz2;

       zt   = [zt1; zt2];
       ZC2  = [ZC2 zt];

       glno=glno+1;
       GL2 = [GL2 glno];

       lz = lz2;

     else  

       xt1  = mesh2d.xc(:,e);
       yt1  = mesh2d.yc(:,e);
   
       xt2  = xt1;
       yt2  = yt1;

       xt = [xt1; xt2];
       yt = [yt1; yt2];

       XC2  = [XC2 xt];
       YC2  = [YC2 yt];

       lz1    = lz + dz;
       lz2    = lz + dz2;
       zt1(1) = lz;
       zt1(2) = lz1;
       zt1(3) = lz1;
       zt1(4) = lz;
  
       zt2(1) = lz2;
       zt2(2) = lz2;
       zt2(3) = lz2;
       zt2(4) = lz2;

       zt   = [zt1; zt2];
       ZC2  = [ZC2 zt];

       glno=glno+1;
       GL2 = [GL2 glno];

       lz=lz2;

     end

     XC2=[];
     YC2=[];
     ZC2=[];
     GL2=[]; 

   end           % k=1:nz

end   % function
%---------------------------------------------------------------------- 

