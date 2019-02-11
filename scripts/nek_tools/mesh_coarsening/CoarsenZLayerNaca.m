function ifcl = CoarsenZLayerNaca(il,nz,nlayers,LX,LY,LE,mesh2d,cz_pl)

% Test 1.
  Zskip = 11;       % Start 'z' refinement after Zskip 2D layer.
 
  nel_lay = length(LE);
  
  ifcl = 0;       % if coarsen layer
  if il==Zskip+1
%   Coarsen entire layer
    ifcl = 1;
  elseif il==Zskip+3
    ifcl = 1;
%  elseif il==Zskip+5
%    ifcl = 1;
  end
  
  if il==nlayers
    ifcl = 0;
  end

% Coarsened layer needs to have multiple of 4 elements for periodicity  
  if mod(nz,4)~=0
    ifcl=0;
  end

% If previous layer was 'z' coarsened  
  if cz_pl
    ifcl = 2;
  end

% We don't modify layers that have been used to modify 'x' resolution
  for j=1:length(LE)
    e=LE(j);
    if strcmpi(mesh2d.EType{e},'e1') || strcmpi(mesh2d.EType{e},'e3') || strcmpi(mesh2d.EType{e},'e4')
      ifcl=0;
    end
  end  


end   %function
%---------------------------------------------------------------------- 

