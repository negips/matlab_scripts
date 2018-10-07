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

% Define an element type  
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
  if nc>0
    disp(['Coarsen ' num2str(nc) ' Elements in layer ' num2str(i)])
  end  

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

% Not done for multiple curvature definitions right now
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
vtkwrite(vfname,'polydata','tetrahedron',xvtk,yvtk,zvtk,polydata,'binary')


% Generate 3D mesh
nz0=4;
Lz=-0.02;
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
vtkwrite(vfname,'polydata','hexahedron',xvtk,yvtk,zvtk,polydata,'binary')
 
glls = [];
for i=1:length(mesh3d.LayerGEl)
  g1 = mesh3d.LayerGEl{i};
  glls = [glls g1];
end
plot(glls)



