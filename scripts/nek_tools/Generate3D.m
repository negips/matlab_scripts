function [NewE3, NewX1, NewY1, NewX2, NewY2, NewBC3, NewCEl3] = Generate3D(NewE,NewX,NewY,NewBC,NewCEl,nlayers,nz0,Lz);

% Generate Layers in 3rd dimension
%---------------------------------------------------------------------- 
%           Element arrangement
%----------------------------------------------------------------------
%
%                         
%                           
%                   x7----------------x6
%                  /|     f6(back)    /|   
%                 /                  / | 
%                /  |               /  |
%               /                  /   |
%     f3(left) /    |  f2(top)    /    |
%             /                  /     |    f1(right) 
%            /     x8- -  - -  -/ - - x5    
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

% Same order for y


% Test 1.
Zskip = 5;       % Start 'z' refinement after Zskip 2D layer.

lez   = NewE;
lxz   = NewX;
lyz   = NewY;
lbcz  = NewBC;
lcelz = NewCEl;

nz = nz0;        % No of elements in 'z' in the first layer.

for il=1:nlayers
 
  LE   = lez{il}; 
  LX   = lxz{il};
  LY   = lyz{il};
  LZ   = zeros(size(lxz{il}));

  LBC  = lbcz{il};  
  LCEl = lcelz{il};

  nel_lay = length(LE);
  maxdlo  = MaxDLO(LX,LY,nel_lay);    
  maxdlv  = MaxDLV(LX,LY,nel_lay);

  ifcl = 0;       % if coarsen layer
  if maxdlo/dz>2 || maxdlv/dz>2 && i>Zskip
%   Coarsen entire layer
    ifcl = 1;
  end 

  if ~ifcl
%   No z coarsening
    dz = Lz/nz;

    for j=1:nz
      for i=1:nel_lay
      
      end
    end       

  end

end


NewE3   = [];
NewX1   = [];
NewX2   = [];
NewY1   = [];
NewY2   = [];
NewBC3  = [];
NewCEl3 = [];


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



