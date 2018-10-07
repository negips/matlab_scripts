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
  Zskip = 5;       % Start 'z' refinement after Zskip 2D layer.
  
  nz = nz0;        % No of elements in 'z' in the first layer.
  dz = Lz/nz;
  
  XC   = [];
  YC   = [];
  ZC   = [];
  GL3D = [];       % Global no for 3d elements
  glno = 0;
  
  il=1;

  cz_pl = 0;       % if previous layer was coarsened
  
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
%     No z coarsening
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
        [XC1,YC1,ZC1,GL1]= QuadExpansion(mesh2d,il,nz,Lz,cz_pl);
  
        il=il+1;
        GL1=GL1+glno;
        glno=glno+length(GL1);
        XC   = [XC XC1];
        YC   = [YC YC1];
        ZC   = [ZC ZC1];
        GL3D = [GL3D GL1];

        if ~isempty(GL1)
          cz_pl=1;
        else
          cz_pl=0;
        end  
  
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




