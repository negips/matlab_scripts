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

  if ifperiodic
    zf1='P  ';
    zf2='P  ';
  else
    disp('Non-periodic BC on the Z faces')
    disp(['Setting zf1=W;    zf2=W'])
    zf1='W  ';
    zf2='W  ';  
  end  


 
  nz = nz0;        % No of elements in 'z' in the first layer.
  dz = Lz/nz;
  
  XC   = [];
  YC   = [];
  ZC   = [];
  GL3D = [];       % Global no for 3d elements
  EF   = [];       % element face bc values

  glno = 0;
  
  il=1;

  cz_pl = 0;       % if previous layer was coarsened
  
  while il<=nlayers
    
    ind  = mesh2d.layerindex{il}; 
    LE   = mesh2d.globalno(ind); 
    LX   = mesh2d.xc(:,ind);
    LY   = mesh2d.yc(:,ind);
    nel_lay = length(LE);
 
    ifcl = CoarsenZLayer(il,nlayers,LX,LY,LE);

    if ~ifcl
%     No z coarsening

      for j=1:nel_lay
        XC1  = [];
        YC1  = [];
        ZC1  = [];
        GL1  = [];       % Global no for 3d elements
        EF1  = [];       % BC on the faces
        e    = LE(j);    % Element number of the 2D mesh
       
        [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_0(e,j,il,nz,Lz,mesh2d,zf1,zf2);

        XC = [XC XC1];
        YC = [YC YC1];
        ZC = [ZC ZC1];
        GL3D = [GL3D GL1];
        EF = [EF EF1];

        LayerGEl{e}=GL1; 
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
        [XC1,YC1,ZC1,GL1,LayerGEl,EF1]= QuadExpansion(mesh2d,LayerGEl,il,nz,Lz,cz_pl,zf1,zf2);
  
        il=il+1;
        GL1=GL1+glno;
        glno=glno+length(GL1);
        XC   = [XC XC1];
        YC   = [YC YC1];
        ZC   = [ZC ZC1];
        GL3D = [GL3D GL1];
        EF   = [EF EF1];

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
  mesh3d.LayerGEl = LayerGEl;
  mesh3d.EF   = EF;

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

function ifcl = CoarsenZLayer(il,nlayers,LX,LY,LE)

% Test 1.
  Zskip = 5;       % Start 'z' refinement after Zskip 2D layer.
 
  nel_lay = length(LE);
  maxdlo  = MaxDLO(LX,LY,nel_lay);    
  maxdlv  = MaxDLV(LX,LY,nel_lay);
  
  ifcl = 0;       % if coarsen layer
  if il==Zskip+1
%   Coarsen entire layer
    ifcl = 1;
%  elseif il==Zskip+4
%    ifcl = 1;
%  elseif il==Zskip+8
%    ifcl = 1;
  end
  
  if il==nlayers
    ifcl = 0;
  end  

end   %function
%---------------------------------------------------------------------- 


