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
Zskip = 2;       % Start 'z' refinement after Zskip 2D layer.

nz = nz0;        % No of elements in 'z' in the first layer.
dz = Lz/nz;

XC = [];
YC = [];
ZC = [];

il=1;

while il<=3 %nlayers
  
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
  if il>Zskip
%   Coarsen entire layer
    ifcl = 1;
  end

  if il==nlayers
    ifcl = 0;
  end  

  [il ifcl]    

  if ~ifcl
%   No z coarsening
    dz = Lz/nz;
    for j=1:nel_lay
      lz = 0;  
      for k=1:nz
%        xt  = [LX(:,j); LX(:,j)];
%        XC  = [XC xt];
%        yt  = [LY(:,j); LY(:,j)];
%        YC  = [YC yt];
%        zt1 = zeros(4,1) + lz;
%        lz  = lz + dz;
%        zt2 = zeros(4,1) + lz;
%        zt  = [zt1; zt2];
%        ZC  = [ZC zt]; 
      end
    end

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
      [XC1,YC1,ZC1,XC2,YC2,ZC2] = QuadExpansion(mesh2d,il,nz,Lz);
      il=il+1;
      XC = [XC XC1];
      YC = [YC YC1];
      ZC = [ZC ZC1];

      continue
    end  

    il=il+1;      % temp
  end % ~ifcl

end   % while il<=nlayers

mesh3d.XC = XC;
mesh3d.YC = YC;
mesh3d.ZC = ZC;


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

function [XC1,YC1,ZC1,XC2,YC2,ZC2] = QuadExpansion(mesh2d,il,nz,Lz);

% Coarsen along z
  il1   = il;  
  ind1  = mesh2d.layerindex{il1}; 
  LE1   = mesh2d.globalno(ind1); 
  LX1   = mesh2d.xc(:,ind1);
  LY1   = mesh2d.yc(:,ind1);
  cbc1  = mesh2d.cbc(:,ind1);

  XC1   = [];
  YC1   = [];
  ZC1   = [];

  XC2   = [];
  YC2   = [];
  ZC2   = [];

  xt1   = zeros(4,1);
  xt2   = zeros(4,1);
  yt1   = zeros(4,1);
  yt2   = zeros(4,1);
  zt1   = zeros(4,1);
  zt2   = zeros(4,1);

  dz = Lz/nz;
  for j=1:6 %nel_lay

    elf1 = cbc1(1,j).connectsto;    % Element no connecting face 1
    elf3 = cbc1(1,j).connectsto;    % Element no connecting face 3
    elf4 = cbc1(4,j).connectsto;    % Element no connecting face 4
%   Element types
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

%   Build the first layer 
    lz = 0;  
    for k=1:nz

      if mod(k,4)==1
        xt  = [LX1(:,j); LX1(:,j)];
        XC1  = [XC1 xt];
        yt  = [LY1(:,j); LY1(:,j)];
        YC1  = [YC1 yt];
        zt1 = zeros(4,1) + lz;
        lz  = lz + dz;
        zt2 = zeros(4,1) + lz;
        zt  = [zt1; zt2];
        ZC1  = [ZC1 zt];

      elseif mod(k,4)==2
        xt1   = LX1(:,j);
        yt1   = LY1(:,j);

        xt2(1) = mesh2d.xc(1,elf4);
        yt2(1) = mesh2d.yc(1,elf4);
        if strcmpi(f1t,'e4')
          xt2(2) = mesh2d.xc(2,elf1);
          yt2(2) = mesh2d.yc(2,elf1);
        else
          xt2(2) = LX1(2,j);
          yt2(2) = LY1(2,j);
        end

        if strcmpi(f3t,'e4')
          xt2(3) = mesh2d.xc(3,elf3);
          yt2(3) = mesh2d.yc(3,elf3);
        else
          xt2(3) = LX1(3,j);
          yt2(3) = LY1(3,j);
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

      elseif mod(k,4)==3

%       face 5            
        xt1(1) = mesh2d.xc(1,elf4);
        yt1(1) = mesh2d.yc(1,elf4);
        if strcmpi(f1t,'e4')
          xt1(2) = mesh2d.xc(2,elf1);
          yt1(2) = mesh2d.yc(2,elf1);
        else
          xt1(2) = LX1(2,j);
          yt1(2) = LY1(2,j);
        end

        if strcmpi(f3t,'e4')
          xt1(3) = mesh2d.xc(3,elf3);
          yt1(3) = mesh2d.yc(3,elf3);
        else
          xt1(3) = LX1(3,j);
          yt1(3) = LY1(3,j);
        end
        xt1(4) = mesh2d.xc(4,elf4);
        yt1(4) = mesh2d.yc(4,elf4);

%       face 6        
        xt2   = LX1(:,j);
        yt2   = LY1(:,j);

        xt    = [xt1; xt2];
        yt    = [yt1; yt2];
        XC1    = [XC1 xt];
        YC1    = [YC1 yt];

        zt1 = zeros(4,1) + lz;
        lz  = lz + dz;
        zt2 = zeros(4,1) + lz;
        zt  = [zt1; zt2];
        ZC1  = [ZC1 zt];

      elseif mod(k,4)==0
        xt  = [LX1(:,j); LX1(:,j)];
        XC1  = [XC1 xt];
        yt  = [LY1(:,j); LY1(:,j)];
        YC1  = [YC1 yt];
        zt1 = zeros(4,1) + lz;
        lz  = lz + dz;
        zt2 = zeros(4,1) + lz;
        zt  = [zt1; zt2];
        ZC1  = [ZC1 zt];

      end         % mod(k,4)

    end           % k=1:nz


%   Build second layer


  end             % nel_layer




end
%---------------------------------------------------------------------- 




