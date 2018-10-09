%     First attempts at modifying the mesh

clear
clc
close all


ifplot = 0;
%casename = 'saab_wing2d';
%casename = 'saab750k';
%casename = 'lu';             % Doesn't work
casename = 'fluent_plus2';
%casename = 'saab10k';
svfname  = [casename '.mat'];

disp(['CaseName: ' casename])

rea = Nek_ReadRea(casename);

n=rea.mesh.nelg;
ndim=rea.mesh.ndim;
cmap = jet(n);

% Find element with bc='v  ' as well as 'O  '
nfaces=2*ndim;
bcels=0;
Elno  = [];
Vface = [];
Oface = [];
OBound = 'O  ';
VBound = 'v  ';
for i=1:n

   ifV=0;
   ifO=0;
   for j=1:nfaces
     bc=rea.mesh.cbc(j,i).bc;
     if strcmpi(bc,VBound)
        ifV=1;
        vface=j;
     end   
     if strcmpi(bc,OBound)
        ifO=1;
        oface=j;
     end
   end   
    
   if (ifV && ifO)
     bcels = bcels+1;
     Elno = [Elno i];
     Vface= [Vface vface];
     Oface= [Oface oface];
   end  

end

% There should only be 2 elements

l1=bcels;
cmap = jet(l1);
%for i=1:bcels
%  glno=Elno(i);
%  xt = rea.mesh.xc(:,glno);
%  yt = rea.mesh.yc(:,glno);  
%  fill(xt,yt,cmap(i,:)); hold on
%end

% Lets take bottom element as starting element
if length(Elno)>1
  gl1 = Elno(1);
  yt1=mean(rea.mesh.yc(:,gl1));
  gl2 = Elno(2);
  yt2=mean(rea.mesh.yc(:,gl2));
  
  ly_el    = [];       % Element number in the layer
  ly_fopO  = [];       % Face no opposite the 'O  ' face
  ly_fopV  = [];       % Face no opposite the 'v  ' face
  
  if yt1<yt2
    ly_el(1) = Elno(1);
    of = Oface(1)+2;
    if of>4
      of=of-4;
    end
    vf = Vface(1)+2;
    if(vf>4)
      vf=vf-4;
    end  
    ly_fopO=of;
    ly_fopV=vf;
  else
    ly_el(1) = Elno(2); 
    of = Oface(2)+2;
    if of>4
      of=of-4;
    end
    vf = Vface(2)+2;
    if(vf>4)
      vf=vf-4;
    end  
    ly_fopO=of;
    ly_fopV=vf;
  end  
elseif (length(Elno)==1)
  gl1 = Elno;
  
  ly_el    = [];       % Element number in the layer
  ly_fopO  = [];       % Face no opposite the 'O  ' face
  ly_fopV  = [];       % Face no opposite the 'v  ' face
  
  ly_el(1) = Elno;
  of = Oface(1)+2;
  if of>4
    of=of-4;
  end
  vf = Vface(1)+2;
  if(vf>4)
    vf=vf-4;
  end  
  ly_fopO=of;
  ly_fopV=vf;

end      

% Build the first layer
done=0;
e=1;
while (~done)
  e1  = ly_el(e);       % current element no
  vf1 = ly_fopV(e);     % face opposite 'v  '
  of1 = ly_fopO(e);     % face opposite 'o  '

  e2  = rea.mesh.cbc(of1,e1).connectsto;        % Next element no.
  of  = rea.mesh.cbc(of1,e1).onface;            % on this face 
  of2 = of+2;
  if (of2>4)
    of2=of2-4;
  end
  vf=0;
  for j=1:nfaces
    bc=rea.mesh.cbc(j,e2).bc;
    if strcmpi(bc,VBound) || strcmpi(bc, 'on ')
      ifV=1;
      vf=j;
    end
  end
  vf2=vf+2;
  if (vf2>4)
    vf2=vf2-4;
  end      

  e=e+1;
  ly_el(e)   = e2;
  ly_fopV(e) = vf2;
  ly_fopO(e) = of2;

  bc = rea.mesh.cbc(of2,e2).bc;
  if ~strcmpi(bc,'E  ')
    done=1;
%   If we don't hit the outflow again, then its not a C-mesh type layer    
    if ~strcmpi(bc,OBound)
      MeshC(1)=0;             % This layer is not Ctype
    else
      MeshC(1)=1;             % This layer is Ctype
    end

  end

end

% Plotting
l1=e;
if ifplot
  cmap = jet(l1);
  for i=1:l1
    glno=ly_el(i);
    xt = rea.mesh.xc(:,glno);
    yt = rea.mesh.yc(:,glno);  
    fill(xt,yt,cmap(i,:)); hold on
  end
end

LayersEl{1}=ly_el;
LayersFopV{1}=ly_fopV;
LayersFopO{1}=ly_fopO;

% Done building first layer

% Build the rest of the layers
finished_layers=0;
nlayers = 1;
while (~finished_layers)
  ply_el   = LayersEl{nlayers};         % Previous Layer's Element number
  ply_fopO = LayersFopO{nlayers};       % Previous Layer's Face opposite the 'O  '
  ply_fopV = LayersFopV{nlayers};       % Previous Layer's Face opposite the 'v  '


  e1  = ply_el(1);      % element no in the previous layer
  vf1 = ply_fopV(1);    % face opposite 'v  '
  e2  = rea.mesh.cbc(vf1,e1).connectsto;
  vf  = rea.mesh.cbc(vf1,e1).onface;
  vf2=vf+2;
  if (vf2>4)
    vf2=vf2-4;
  end  

  of=0;
  for j=1:nfaces
    bc=rea.mesh.cbc(j,e2).bc;
    if strcmpi(bc,OBound)
      of=j;
    end
  end
  of2=of+2;
  if (of2>4)
    of2=of2-4;
  end  

  ly_el    = e2;
  ly_fopO  = of2;
  ly_fopV  = vf2;

  done=0;

  e=1;
  while (~done)
    e1  = ly_el(e);
    of1 = ly_fopO(e);
    vf1 = ly_fopV(e);  

    e2  = rea.mesh.cbc(of1,e1).connectsto;        % Next element no.
    of  = rea.mesh.cbc(of1,e1).onface;            % on this face 
    of2 = of+2;
    if (of2>4)
      of2=of2-4;
    end

    plel   = ply_el(e+1);
    plfopv = ply_fopV(e+1); 
    e3 = rea.mesh.cbc(plfopv,plel).connectsto;
    if e3 ~= e2
      disp(['In consistent elements: ', num2str(nlayers), '; e: ',num2str(e)])
      disp(['e2:',num2str(e2), ';  e3:',num2str(e3)])
      vf2=0;
    else
      vf = rea.mesh.cbc(plfopv,plel).onface;
      vf2=vf+2;
      if vf2>4
        vf2=vf2-4;
      end
    end  
  
    e=e+1;
    ly_el(e)   = e2;
    ly_fopV(e) = vf2;
    ly_fopO(e) = of2;
  
    bc = rea.mesh.cbc(of2,e2).bc;
    if ~strcmpi(bc,'E  ')
      done=1;
    end
  
  end       % ~done
 
  nlayers = nlayers+1;
  LayersEl{nlayers}   = ly_el;         % New Layer's Element numbers
  LayersFopO{nlayers} = ly_fopO;       % New Layer's Face opposite the 'O  '
  LayersFopV{nlayers} = ly_fopV;       % New Layer's Face opposite the 'v  '

% If we don't hit the outflow again, then its not a C-mesh type layer    
  if ~strcmpi(bc,OBound)
    MeshC(nlayers)=0;             % This layer is not Ctype
  else
    MeshC(nlayers)=1;             % This layer is Ctype
  end

% Test for end of layers
  el2=ly_el(1);
  fv2=ly_fopV(1);
  bc2=rea.mesh.cbc(fv2,el2).bc;

  if ~strcmpi(bc2,'E  ')
    finished_layers=1;
  else
    el3=rea.mesh.cbc(fv2,el2).connectsto;
    % Check all layers
    for kk=nlayers-1:-1:1
      tmp = find(LayersEl{kk}==el3);
      if ~isempty(tmp)
        finished_layers=1;
        break
      end
    end  
  end  



  % Plotting
  l1=e;
  if ifplot
    cmap = jet(l1);
    for i=1:l1
      glno=ly_el(i);
      xt = rea.mesh.xc(:,glno);
      yt = rea.mesh.yc(:,glno);  
      fill(xt,yt,cmap(i,:)); hold on
    end
  end

end   % ~finished_layers 

%arrange_matrices

[LayerX,LayerY,LayerE,LayerBC,LayerCEl] = ArrangeMatrices(nlayers,LayersEl,LayersFopO,LayersFopV,MeshC,rea,ndim);

clearvars -except nlayers LayersEl LayersFopO LayersFopV LayerX LayerY LayerE LayerBC LayerCEl MeshC rea n ndim svfname
save(svfname)



