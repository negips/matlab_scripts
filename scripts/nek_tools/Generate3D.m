function [mesh3d] = Generate3D(mesh2d,nlayers,nz0,Lz,ifperiodic);

% Generate Layers in 3rd dimension
%---------------------------------------------------------------------- 
%           Element arrangement
%----------------------------------------------------------------------
%
%                         
%                           
%                   x7-----------------x6
%                  /|     f6(back)    /|   
%                 /                  / | 
%                /  |               /  |
%               /                  /   |
%     f3(left) /    |  f2(top)    /    |
%             /                  /     |    f1(right) 
%            /      x8- -  - - -/ - - -x5    
%           x3-----------------x2     /  
%           |     /            |     /
%           |        f4(bottom)|    / 
%           |   /              |   /
%           |                  |  /
%           | /                | /  
%           x4-----------------x1
%                 
%                 f5(front) 
%
%
%---------------------------------------------------------------------- 
%---------------------------------------------------------------------- 

% Initial 2d Layer:
% (x1, x2, x3, x4)
% (y1, y1, y3, y4)
% Faces: (f1, f2, f3, f4)


% Defined 3D face coordinates:

% f1: (x1, x2, x6, x5)
% f2: (x2, x3, x7, x6)
% f3: (x3, x4, x7, x8)
% f4: (x4, x1, x5, x8)
% f5: (x1, x2, x3, x4)
% f6: (x5, x6, x7, x8)

% Same order for y,z


% Test 1.
Zskip = 10;       % Start 'z' refinement after Zskip 2D layer.

nz = nz0;        % No of elements in 'z' in the first layer.
dz = Lz/nz;

XC   = [];
YC   = [];
ZC   = [];
GL3D = [];       % Global no for 3d elements
glno = 0;

il=1;

while il<=nlayers
  
  ind  = mesh2d.layerindex{il}; 
  LE   = mesh2d.globalno(ind); 
  LX   = mesh2d.xc(:,ind);
  LY   = mesh2d.yc(:,ind);
  LZ   = zeros(size(LX));
  cbc  = mesh2d.cbc(:,ind);

  nel_lay = length(LE);
  maxdlo  = MaxDLO(LX,LY,nel_lay);    
  maxdlv  = MaxDLV(LX,LY,nel_lay);

  ifcl = 0;       % if coarsen layer
%  if maxdlo/dz>2 || maxdlv/dz>2 && il>Zskip
  if il==Zskip+1
%   Coarsen entire layer
    ifcl = 1;
  end

  if il==nlayers
    ifcl = 0;
  end  

  if ~ifcl
%   No z coarsening
%    if il==5 || il==6 
      dz = Lz/nz;
      for j=1:nel_lay
        lz = 0;  
        for k=1:nz
          xt  = [LX(:,j); LX(:,j)];
          XC  = [XC xt];
          yt  = [LY(:,j); LY(:,j)];
          YC  = [YC yt];
          zt1 = zeros(4,1) + lz;
          lz  = lz + dz;
          zt2 = zeros(4,1) + lz;
          zt  = [zt1; zt2];
          ZC  = [ZC zt];
          glno = glno+1;
          GL3D = [GL3D glno];
        end       % k
      end         % j
%    end           % il  
    il=il+1;

  else

%   Coarsen along z
    ifodd = mod(nz,2);
    ifquad = 0;
    if ~ifodd
      rem = mod(nz,4);
      if rem==0
        ifquad = 1;
      else
        ifquad = 0;
      end
    end

    if ifquad
      [XC1,XC2,XC3,XC4,YC1,YC2,YC3,YC4,ZC1,ZC2,ZC3,ZC4,GL1,GL2,GL3,GL4] = QuadExpansion(mesh2d,il,nz,Lz);

      il=il+2;
      GL1=GL1+glno;
      GL2=GL2+glno;
      glno=glno+length(GL1)+length(GL2);
      XC   = [XC XC1 XC2];
      YC   = [YC YC1 YC2];
      ZC   = [ZC ZC1 ZC2];
      GL3D = [GL3D GL1 GL2];

      nz = nz/2;

      continue
    end 

  end % ~ifcl

end   % while il<=nlayers

mesh3d.XC   = XC;
mesh3d.YC   = YC;
mesh3d.ZC   = ZC;
mesh3d.GL3D = GL3D;


end   % function

%---------------------------------------------------------------------- 

function maxdlo = MaxDLO(LX,LY,nel)

%   Length of sides 1 and 3   % Facing 'O'
    dlo = 0;  
    for j=1:nel  
      l1o = sqrt( (LX(1,j)-LX(2,j))^2 + (LY(1,j) - LY(2,j))^2);
      l2o = sqrt( (LX(3,j)-LX(4,j))^2 + (LY(3,j) - LY(4,j))^2);
      dlo = max([l1o l2o dlo]);
    end
    maxdlo = dlo;  

end % function
%---------------------------------------------------------------------- 

function maxdlv = MaxDLV(LX,LY,nel)
          
%   Length of sides 2 and 4   % Facing 'V'
    dlv=0;  
    for j=1:nel  
      l1v = sqrt( (LX(2,j)-LX(3,j))^2 + (LY(2,j) - LY(3,j))^2);
      l2v = sqrt( (LX(4,j)-LX(1,j))^2 + (LY(4,j) - LY(1,j))^2);
      dlv = max([l1v l2v dlv]);
    end
    maxdlv = dlv;


end % function
%---------------------------------------------------------------------- 

function [XC1,XC2,XC3,XC4,YC1,YC2,YC3,YC4,ZC1,ZC2,ZC3,ZC4,GL1,GL2,GL3,GL4] = QuadExpansion(mesh2d,il,nz0,Lz);

% Coarsen along z
  il1   = il;  
  ind1  = mesh2d.layerindex{il1}; 
  LE1   = mesh2d.globalno(ind1); 

  nel_lay = length(ind1);

  XC1 = []; YC1 = []; ZC1 = []; GL1 = [];

  XC2 = []; YC2 = []; ZC2 = []; GL2 = [];

  XC3 = []; YC3 = []; ZC3 = []; GL3 = [];

  XC4 = []; YC4 = []; ZC4 = []; GL4 = [];
 
  j=1;
  while j<=nel_lay

    e=LE1(j);
    nz=nz0;

    etype = mesh2d.EType(e);

    if strcmpi(etype,'s') 
      [XC1,XC2,XC3,XC4,YC1,YC2,YC3,YC4,ZC1,ZC2,ZC3,ZC4,GL1,GL2,GL3,GL4,j]=BuildLayer_S(XC1,XC2,XC3,XC4,YC1,YC2,YC3,YC4,ZC1,ZC2,ZC3,ZC4,GL1,GL2,GL3,GL4,e,j,il,nz,Lz,mesh2d);
    end

  end             % nel_layer


end
%---------------------------------------------------------------------- 

function [XC1,XC2,XC3,XC4,YC1,YC2,YC3,YC4,ZC1,ZC2,ZC3,ZC4,GL1,GL2,GL3,GL4,j]=BuildLayer_S(XC1,XC2,XC3,XC4,YC1,YC2,YC3,YC4,ZC1,ZC2,ZC3,ZC4,GL1,GL2,GL3,GL4,e,j,il,nz,Lz,mesh2d);

%  Build the first layer
   xt1 = zeros(4,1); yt1 = zeros(4,1); zt1 = zeros(4,1);
   xt2 = zeros(4,1); yt2 = zeros(4,1); zt2 = zeros(4,1);
   xt3 = zeros(4,1); yt3 = zeros(4,1); zt3 = zeros(4,1);
   xt4 = zeros(4,1); yt4 = zeros(4,1); zt4 = zeros(4,1);

   elf1 = mesh2d.cbc(1,e).connectsto;    % Element no connecting face 1
   elf3 = mesh2d.cbc(1,e).connectsto;    % Element no connecting face 3
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
         xt2(3) = mesh2d.xc(3,j);
         yt2(3) = mesh2d.yc(3,j);
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

%      face 5            
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

%      face 6        
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

   end           % k=1:nz


end
%---------------------------------------------------------------------- 






