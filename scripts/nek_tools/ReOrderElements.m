function [Mesh2D] = ReOrderElements(NewE,NewX,NewY,NewBC,NewCEl,oldmesh,cdef);

nlayers = length(NewE);
tnel=0;
NewGN = [];
OldGN = [];
XC    = [];
YC    = [];
glno=0;
LayerEnd = [];
for il=1:nlayers
  EL = NewE{il};
  X  = NewX{il};
  Y  = NewY{il};
  nel = length(EL);
  tnel = tnel+nel;
  for j=1:nel
    glno = int32(glno+1);
    OldGN = [OldGN EL(j)];
    NewGN = [NewGN glno];
    XC    = [XC X(:,j)];
    YC    = [YC Y(:,j)];
  end
  LayerEnd = [LayerEnd glno]; 

end

Mesh2D.nelg=tnel;
Mesh2D.globalno = NewGN;
groupno = int32(zeros(1,tnel));
Mesh2D.groupno = groupno;
Mesh2D.xc = XC;
Mesh2D.yc = YC;

% Assuming I don't touch the curved faces
Mesh2D.ncurve=oldmesh.ncurve;
[C,IA,IB]=intersect(oldmesh.curveieg,OldGN,'stable');
ncurve=length(IB);

curveface = [];
for ic=1:ncurve
  curveieg(ic)=NewGN(IB(ic));
  curveparams(:,ic) = oldmesh.curveparams(:,ic);

% Find the face number for the curve
  nfaces = 4;
  ieg = curveieg(ic);
  iffound = 0;
  Layerfound = 0;
  ieg0=0;
  for i=1:nlayers
    if ieg<=LayerEnd(i)
%     Element belongs to this layer
      j=ieg-ieg0;
      Layerfound=i;
      for k=1:nfaces
        bc = NewBC{i}{j}(k,:);
        if strcmpi(bc,cdef)
          curveface = [curveface k];
          iffound = 1;    
          break
        end
        if iffound
          break
        end
      end  % k=1:nfaces
    end    % ieg<=LayerEnd
    ieg0=LayerEnd(i);

    if (iffound)
      break
    end  
  end      % i=1:nlayers
  if ~iffound
    disp(['Element with No ',num2str(ieg) ' and definition ', cdef ,' not found in Layer ', num2str(Layerfound)])
  end  

end
Mesh2D.curveieg=curveieg;
Mesh2D.curveparams=curveparams;
Mesh2D.curveface=curveface;


% Make cbc structure
nfaces=4;
for il=1:nlayers
  nel = length(NewE{il});
  for j=1:nel
    for k=1:nfaces
      bc = NewBC{il}{1}(k,:);



  end

end   % il=1:nlayers

end   % end function


