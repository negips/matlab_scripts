% First attempts at coarsening the mesh

clear
clc
close all


%           Element arrangement
%
%
%                 f2
%           x3-----------x2      
%           |            |
%           |            |
%         f3|            |f1 'O  '
%           |            |
%           |            |
%           x4-----------x1
%                 f4
%               'v  '


%load saab_wing2d.mat
%load saab750k.mat
load fluent_plus2.mat

skiplayers = 2;         % Need to skip first layer since its smaller than the others
ifplot = 0;
if3dplot = 0;

disp(['Total Number of Elements: ', num2str(rea.mesh.nelg)])

for i=1:nlayers
  j=nlayers-i+1;
  OldE{i}=LayerE{j};
  OldX{i}=LayerX{j};
  OldY{i}=LayerY{j};
  OldBC{i}=LayerBC{j};
  OldCEl{i}=LayerCEl{j};

  NewE{i}=OldE{i};
  NewX{i}=OldX{i};
  NewY{i}=OldY{i};
  NewBC{i}=OldBC{i};
  NewCEl{i}=OldCEl{i};
  NewCoF{i}=ConnectedOnFace(NewCEl{i});

  l1=length(NewE{i});
  for j=1:l1
    NewET{i}{j}='s';
  end  

end


% Test coarsening algorithm 1
% Coarsen Layer by Layer
% Define aspect ratio as 'O' face lengths to 'V' face lengths
layer_start = skiplayers+1;
new_nelg = 0;
if ifplot
  fig2=figure(2);
else
  fig2=[];
end  

for i=1:nlayers

  LE=NewE{i};
  LX=NewX{i};
  LY=NewY{i};
  LCEl=NewCEl{i}; 

  if (i<layer_start)
    nels_layer = length(LX);    
    new_nelg = new_nelg+nels_layer; 

    if (ifplot)
      cmap = jet(nels_layer);
      for j=1:nels_layer 
        xt = LX(:,j);
        yt = LY(:,j);

        figure(2)
        fill(xt,yt,cmap(j,:)); hold on
      end 
    end

    continue
  end  

  l1 =length(LX);
  if i==layer_start
    iflocked=zeros(l1,1);
  end  
  cmap = jet(l1);
  l2 = length(OldX{i});
  ifc = zeros(l2,1); 
  for j=1:l1
      
%    ifc(j) = CoarsenCriteria(LX,LY,j,i,iflocked);
    ifc(j) = CoarsenKDJ(LX,LY,j,i,iflocked);
%    ifc(j) = CoarsenSaab(LX,LY,j,i,iflocked);
      
    xt = LX(:,j);
    yt = LY(:,j);

%    figure(1)
%    fill(xt,yt,cmap(j,:)); hold on
  end

  if i==nlayers
    ifc(:)=0;
  end  

  c_ind = find(ifc);
  nc = length(c_ind);
  disp(['Coarsen ' num2str(nc) ' Elements in layer ' num2str(i)])

  ifplot =0;
% Coarsen layer in consecutive pairs
  [LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,NewET] = Coarsen2DLayer(LE,LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,NewET,i,fig2,ifplot);
% Modify all subsequent Layers
  if3skip=1;
  [NewE, NewX, NewY, NewBC, NewCEl, NewCoF, NewET, iflocked]=Create2DLayers(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,NewET,i,ifc,if3skip,fig2,ifplot);
  nels_layer = length(NewX{i});
  new_nelg = new_nelg+nels_layer;
end
disp(['Total Number of Elements: ', num2str(new_nelg)])

% Not done for multiple definitions right now
curvedef= 'W  ';
mesh2d = ReOrderElements(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,NewET,rea.mesh,curvedef); 

polydata = [];

nel=mesh2d.nelg;
for i=1:nel

  xgll(:,i) = mesh2d.xc(:,i);
  ygll(:,i) = mesh2d.yc(:,i);

  p0 = (2^2)*(i-1);

  f1 = [0 1 2 3] + 1 + p0;

  polydata = [polydata; f1];

end  

xvtk=xgll(:);
yvtk=ygll(:);
zvtk=0*xvtk;

vfname = 'test2d.vtk';
vtkwrite(vfname,'polydata','tetrahedron',xvtk,yvtk,zvtk,polydata)



% Generate 3D mesh
nz0=16;
Lz=0.02;
ifperiodic=1;
[mesh3d] = Generate3D(mesh2d,nlayers,nz0,Lz,ifperiodic);
[~, nel]=size(mesh3d.XC)

if if3dplot
  fig3 = figure(3); hold on
  for i=1:nel
    Plot3DElement(mesh3d,i,fig3);
  end
end  


% Output as a Nek field file
ndim=3;
N=1;
nel=nel;
nps=0;
ifx=1;
ifu=1;
ifp=0;
Glno=mesh3d.GL3D;
U = [];
V = [];
W = [];
P = []; 
T = []; %mesh3d.GL3D;

P = zeros(2^ndim,nel);
U = zeros(2^ndim,nel);
V = zeros(2^ndim,nel);
W = zeros(2^ndim,nel);
xgll = zeros(2^ndim,nel);
ygll = zeros(2^ndim,nel);
zgll = zeros(2^ndim,nel);

polydata = [];

for i=1:nel

  xgll(:,i) = mesh3d.XC(:,i);
  ygll(:,i) = mesh3d.YC(:,i);
  zgll(:,i) = mesh3d.ZC(:,i);

  p0 = (2^ndim)*(i-1);

  f1 = [0 1 5 4] + 1 + p0;
  f2 = [1 5 6 2] + 1 + p0;
  f3 = [2 6 7 3] + 1 + p0;
  f4 = [3 7 4 0] + 1 + p0;
  f5 = [0 1 2 3] + 1 + p0;
  f6 = [4 5 6 7] + 1 + p0;

  polydata = [polydata; f1; f2; f3; f4; f5; f6];

end  

xvtk=xgll(:);
yvtk=ygll(:);
zvtk=zgll(:);

vfname = 'test.vtk';
vtkwrite(vfname,'polydata','hexahedron',xvtk,yvtk,zvtk,polydata)
 
%---------------------------------------------------------------------- 
function Plot3DElement(mesh3d,i,fig)

  figure(fig)
  
  ind=[1 2 3 4 1];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  ind=[5 6 7 8 5];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  ind=[1 5 6 2 1];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  ind=[4 8 7 3 4];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  xlabel('x')
  ylabel('y')
  zlabel('z')
  view(3)

end
%---------------------------------------------------------------------- 

function newcof = ConnectedOnFace(NewCEl)

  [r nel] = size(NewCEl);
 
  newcof = []; 
  nfaces=r;
  for i=1:nel
    f=int32(zeros(nfaces,1));
    for j=1:nfaces
      if NewCEl(j,i)==0   % Boundary element
        f(j)=0;
      else
        f2=j+2;
        if f2>4
          f2=f2-4;
        end
        f(j)=f2;
      end     % if
  
    end       % j=1:nfaces
    newcof = [newcof f]; 
  
  end         % i=1:nel

end   % function  
%---------------------------------------------------------------------- 

function ifc = CoarsenCriteria(LX,LY,j,i,iflocked)

    ARcut = 2.5;             % Coarsen if Aspect ratio is larger than this.
    start_layer = 5;

%   Find aspect ratio of elements 

%   Length of sides 1 and 3   % Facing 'O'
    l1o = sqrt( (LX(1,j)-LX(2,j))^2 + (LY(1,j) - LY(2,j))^2);
    l2o = sqrt( (LX(3,j)-LX(4,j))^2 + (LY(3,j) - LY(4,j))^2);
    dlo = mean([l1o l2o]);
     
%   Length of sides 2 and 4   % Facing 'V'
    l1v = sqrt( (LX(2,j)-LX(3,j))^2 + (LY(2,j) - LY(3,j))^2);
    l2v = sqrt( (LX(4,j)-LX(1,j))^2 + (LY(4,j) - LY(1,j))^2);
    dlv = mean([l1v l2v]);

    l_ar = dlo/dlv;
    lmax = max([l1o l2o l1v l2v]);

    ifc=0;
    if l_ar>ARcut
      ifc=1; 
    end

    xmid = mean(LX(:,j));
    ymid = mean(LY(:,j));
    rad = sqrt(xmid^2 + ymid^2);  

    if (xmid<0.0 && rad>0.1 )
      ifc=0;
    end

    if xmid>1 && l_ar>1.25
      ifc=1;
    end

%   Maximum length  
    if dlv>0.1
      ifc=0;
    end

%   For the radially emerging elements I refine by number of layers
    if i==start_layer
      theta=atan(ymid/(xmid-0.25))*180/pi;
      if xmid<0.25 && abs(theta)<15
        ifc=1;
      end
    elseif i<=start_layer+4
      theta=atan(ymid/(xmid-0.25))*180/pi;
      if xmid<0.25 && abs(theta)>15 && abs(theta)<75
        ifc=1;
      end
    else
      if xmid<0.25 
        ifc=0;
      end
    end
    
%   Skip first n layers      
    if i<start_layer
      ifc=0;
    end  

    [pts nels]=size(LX);  
%   End condition           
    if (iflocked(j) || j==1 || j>=nels-1)
      ifc=0;
    end
    
%   If this layer has anything locked, lock the whole layer
    if max(iflocked)
      ifc=0;
    end  


end % function
%---------------------------------------------------------------------- 






