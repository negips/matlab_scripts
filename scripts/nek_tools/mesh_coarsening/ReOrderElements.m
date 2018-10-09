function [Mesh2D] = ReOrderElements(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,NewET,oldmesh,cdef);

  nlayers = length(NewE);
  tnel=0;
  NewGN = [];
  OldGN = [];
  XC    = [];
  YC    = [];
  EType = [];
  glno=0;
  LayerEnd = [];
  for il=1:nlayers
    EL = NewE{il};
    X  = NewX{il};
    Y  = NewY{il};
    ET = NewET{il};
    nel = length(EL);
    tnel = tnel+nel;
    lindex = [];
    for j=1:nel
      glno   = int32(glno+1);
      lindex = [lindex glno];
      OldGN  = [OldGN EL(j)];
      NewGN  = [NewGN glno];
      XC     = [XC X(:,j)];
      YC     = [YC Y(:,j)];
      EType  = [EType ET(j)];
    end
    LayerEnd = [LayerEnd glno];
    LayerIndex{il} = lindex; 
  end

  [oldsort ind]=sort(OldGN);
  newsort = NewGN(ind);
  j=1;
  for i=1:oldmesh.nelg
    if i==oldsort(j)
       Old2New(i) = newsort(j);
       j=j+1;
    else
       Old2New(i) = -1;
    end
  end
   
  Mesh2D.nelg=tnel;
  Mesh2D.ndim=oldmesh.ndim;
  Mesh2D.globalno = NewGN;
  groupno = int32(zeros(1,tnel));
  Mesh2D.groupno = groupno;
  Mesh2D.xc = XC;
  Mesh2D.yc = YC;
  Mesh2D.EType = EType;
  
% Assuming I don't touch the curved faces
  [curveieg,curveface,curveparams,ctype,ncurve] = Build2DCurved(oldmesh,Old2New,nlayers,NewBC,LayerEnd,cdef);

  Mesh2D.ncurve=ncurve;
  Mesh2D.curveface=curveface;
  Mesh2D.curveieg=curveieg;
  Mesh2D.curveparams=curveparams;
  Mesh2D.curvetype=ctype;

  
% Make cbc structure
  cbc = Build2DCBC(NewBC,NewCEl,NewCoF,Old2New,NewE,nlayers);
 
  Mesh2D.cbc=cbc; 
  Mesh2D.xfac=oldmesh.xfac;
  Mesh2D.yfac=oldmesh.yfac;
  Mesh2D.xzero=oldmesh.xzero;
  Mesh2D.yzero=oldmesh.yzero;
  Mesh2D.ifre2=oldmesh.ifre2;
  Mesh2D.ifgtp=oldmesh.ifgtp;
  Mesh2D.layerindex=LayerIndex;

end   % end function

%----------------------------------------------------------------------

function [curveieg,curveface,curveparams,ctype,ncurve] = Build2DCurved(oldmesh,Old2New,nlayers,NewBC,LayerEnd,cdef)

  newcurveieg = Old2New(oldmesh.curveieg);
      
  ncold=length(newcurveieg);
 
  curveieg    = []; 
  curveface   = [];
  curveparams = [];
  ctype       = [];
  ncurve = 0;
  for ic=1:ncold
  
%   Find the edge number for the curve
    nfaces = 4;
    newieg = newcurveieg(ic);
    if newieg<0
      disp('Old curved element has been removed. Removing from curved sides')
      continue
    end    

    iffound = 0;
    Layerfound = 0;
    for i=1:nlayers
      if newieg<=LayerEnd(i)
        Layerfound=i;
        break    
      end
    end  

%   Element belongs to this layer
    if Layerfound==0
      disp('Something went wrong in finding new curved element.')
      continue
    elseif Layerfound==1
      i=1;
      ieg0=0;      
    else
      i=Layerfound;
      ieg0=LayerEnd(Layerfound-1);
    end
          
    j=newieg-ieg0;
    for k=1:nfaces
      bc = NewBC{i}{j}(k,:);
      if strcmpi(bc,cdef)
        curveieg  = [curveieg newieg];    
        curveface = [curveface k];
        iffound = 1;
        ncurve=ncurve+1;
        curveieg(ncurve)=newieg;
        curveparams(:,ncurve) = oldmesh.curveparams(:,ic);
        ctype{ncurve}=oldmesh.curvetype{ic};
  
        break
      end
    end  

    if ~iffound
      ed = [];
      for k=1:4
        ed = [ed NewBC{i}{j}(k,:)];
      end    
      disp(['Definition ', cdef, 'not found. Element ', num2str(newieg),' Edges: ',ed])
    end
  
  end
 
end   % function
%----------------------------------------------------------------------
function cbc = Build2DCBC(NewBC,NewCEl,NewCoF,Old2New,NewE,nlayers)

  cbc = [];
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
%        connects = [connects c2el];
        if c2el==0
          connects=0;
        else  
          connects = Old2New(c2el);
        end

        cbc(k,glno).bc = bc;
        cbc(k,glno).connectsto = connects;

        onface = NewCoF{il}(k,j);
        cbc(k,glno).onface = onface;

%       Hard coded since I have not saved information on this        
        cbc(k,glno).param1 = 0;
        cbc(k,glno).param2 = 0;
        cbc(k,glno).param3 = 0;

        ifE  = strcmpi(bc, 'E  ');
        if ~ifE && c2el~=0
          disp(['Inconsistency between BC and connectsto arrays on (il,j,k)=(',num2str(il),num2str(j),num2str(k),')']); 
        end
  
      end
  
    end
  
  end   % il=1:nlayers

end   % function
%---------------------------------------------------------------------- 
