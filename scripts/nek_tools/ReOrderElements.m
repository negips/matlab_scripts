function [Mesh2D] = ReOrderElements(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,oldmesh,cdef);

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
  Mesh2D.ndim=oldmesh.ndim;
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
  ncurve2 = 0;
  for ic=1:ncurve
  
  %  curveieg(ic)=NewGN(IB(ic));
  %  curveparams(:,ic) = oldmesh.curveparams(:,ic);
  
  % Find the face number for the curve
    nfaces = 4;
    ieg = NewGN(IB(ic));    
  %  ieg = curveieg(ic);
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
            ncurve2=ncurve2+1;
            curveieg(ncurve2)=ieg;
            curveparams(:,ncurve2) = oldmesh.curveparams(:,ic);
            ctype{ncurve2}=oldmesh.curvetype{ic};
  
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
  
  Mesh2D.ncurve=ncurve2;
  Mesh2D.curveface=curveface;
  Mesh2D.curveieg=curveieg;
  Mesh2D.curveparams=curveparams;
  Mesh2D.curvetype=ctype;
  
  
  % Make cbc structure
  nfaces=4;
  connects = [];
  glno=0; 
  for il=1:nlayers
    nel = length(NewE{il});
    for j=1:nel
      glno = glno+1;  
      for k=1:nfaces
        bc = NewBC{il}{j}(k,:);
        c2el = NewCEl{il}(k,j);
        connects = [connects c2el];
        ifE  = strcmpi(bc, 'E  ');
        if ~ifE && c2el~=0
          disp(['Inconsistency between BC and connectsto arrays on (il,j,k)=(',num2str(il),num2str(j),num2str(k),')']); 
        end
  
      end
  
    end
  
  end   % il=1:nlayers

  [C,IA,IB]=intersect(connects,OldGN,'stable');
  Cn2 = NewGN(IB);

  s=0;
  connects2 = connects;
  for i=1:length(OldGN)
    ind=find(connects==OldGN(i));
    s=s+length(ind);
    connects2(ind)=NewGN(i);
  end

  ind0 = find(connects==0);
  s0   = length(ind0);
  connects2(ind0)=0;

  s1=length(connects2);
  if s+s0~=length(connects2)
    disp(['Some connecting element numbers are missing ', num2str(s), ' ', num2str(s1), ' ', num2str(s0)])
  end  


  glno=0;
  faceno=0;
  for il=1:nlayers
    nel = length(NewE{il});
    for j=1:nel
      glno = glno+1; 
      for k=1:nfaces
        bc = NewBC{il}{j}(k,:);
        cbc(k,glno).bc = bc;

        faceno=faceno+1;
        cbc(k,glno).connectsto = connects2(faceno);

        onface = NewCoF{il}(k,j);
        cbc(k,glno).onface = onface;

%       Hard coded since I have not saved information on this        
        cbc(k,glno).param1 = 0;
        cbc(k,glno).param2 = 0;
        cbc(k,glno).param3 = 0;

      end
    end
  end   % il=1:nlayers

  Mesh2D.cbc=cbc; 
  Mesh2D.xfac=oldmesh.xfac;
  Mesh2D.yfac=oldmesh.yfac;
  Mesh2D.xzero=oldmesh.xzero;
  Mesh2D.yzero=oldmesh.yzero;
  Mesh2D.ifre2=oldmesh.ifre2;
  Mesh2D.ifgtp=oldmesh.ifgtp;

end   % end function


