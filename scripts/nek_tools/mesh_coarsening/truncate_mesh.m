%     First attempts at modifying the mesh

clear
clc
close all
%
%ifplot = 0;
casename = 'saab750k';
%
modify_mesh2

% Identify global numbers to keep

Keep = LayerE;
for i=1:nlayers
  Keep{i}(:)=0;
end  


% Find Element Numbers to keep
% On the Lower side, find element no closest to x=0.2
stl=38;                   % First full layer
i=stl;
nel = length(LayerE{i});
j=1;
xmid = mean(LayerX{i}(:,j));
while xmid>0.2 && j<nel
  j=j+1;
  xmid = mean(LayerX{i}(:,j));
end

Ls = j;    % This is the starting no to keep for each layer
while xmid<0.8 && j<nel
  Keep{i}(j)=1;
  j=j+1;
  xmid = mean(LayerX{i}(:,j));
end
Le = j-1;    % This is the last no to keep for each layer
enl = 10;    % end layer
for i=stl-1:-1:enl
  Keep{i}(Ls:Le)=1;
end 


fig10=figure(1); hold on
% Old and New Global numbers of the elements
GLOld=rea.mesh.globalno;
GLNew=GLOld;
gl=0;
gg=0;
for i=nlayers:-1:1
  for j=1:length(Keep{i})
    gg=gg+1;
    GLOld(gg)=LayerE{i}(j);
    if Keep{i}(j)>0
      gl = gl+1;
      GLNew(gg) = gl;
%      Plot2DElement(rea.mesh,gg,fig10)
    else
      GLNew(gg) = 0;
      if j==Ls-1 || j==Le+1
        GLNew(gg) = -1;       % Elements connecting to these specified as outflow
      end
      if i==enl-1 && (j>=Ls && j<= Le) 
        GLNew(gg) = -2;       % Elements connecting to these specified as Dirichlet
      end
    end       % Keep>0 
  end
end  
NelgNew = gl;


% I don't want to change orientations.
% So I will modify the old rea.mesh structure

[GLOld,I] = sort(GLOld);
GLNew = GLNew(I);

%[glnew_sorted,I2] = sort(GLNew);

% Modify CBC
nelgold = rea.mesh.nelg;
nfaces  = 2*ndim;

cbco = rea.mesh.cbc;

cbc = rea.mesh.cbc;
outflow = 0;
dirich  = 0;
c2max = 0;
for i=1:nelgold
% If we keep this element      
  if GLNew(i)>0
    for j=1:nfaces
      cto = cbc(j,i).connectsto;
      if cto>0      
        c2new = GLNew(cto);
        c2max = max(c2new,c2max);
      else
        c2new = 0;
      end        
            
      if c2new==-1    % outflow
        outflow=outflow+1;
        cbc(j,i).connectsto=0;
        cbc(j,i).onface=0;
        cbc(j,i).bc='O  ';     % set bc as outflow
%        disp([num2str(outflow), ' Switched to outflow ' num2str(i), ' ', num2str(GLNew(i))])
      elseif c2new==-2      % Dirichlet
        dirich = dirich + 1;
        cbc(j,i).connectsto=0;
        cbc(j,i).onface=0;
        cbc(j,i).bc='v  ';     % set bc as outflow
%        disp([num2str(dirich), ' Switched to Dirichlet ' num2str(i), ' ', num2str(GLNew(i))])
      else
        cbc(j,i).connectsto=c2new;
      end          
    end     % j=1:nfaces
  end       % if (GLNew(i)>0)
end         % i=1:nelgold
disp([num2str(outflow), ' Elements switched to outflow'])
disp([num2str(dirich),  ' Elements switched to Dirichlet'])

% Change Curved Element numbers
crv_ind = [];
for i=1:rea.mesh.Ncurve
  elno = rea.mesh.curveieg(i);
  if GLNew(elno)>0
    crv_ind = [crv_ind i];
    rea.mesh.curveieg(i)=GLNew(i);
  end
end  
crv_ind = sort(crv_ind);

% Throw away unwanted elements
IndF = GLNew<=0;
IndT = find(GLNew>0);
GLNew2= GLNew(IndT);

rea.mesh.xc  = rea.mesh.xc(:,IndT);
rea.mesh.yc  = rea.mesh.yc(:,IndT);
rea.mesh.cbc = cbc(:,IndT);
rea.mesh.curveedge   = rea.mesh.curveedge(crv_ind);
rea.mesh.curveieg    = rea.mesh.curveieg(crv_ind);
rea.mesh.curveparams = rea.mesh.curveparams(:,crv_ind);
rea.mesh.curvetype   = rea.mesh.curvetype(crv_ind);
rea.mesh.Ncurve      = length(crv_ind);
rea.mesh.globalno    = GLNew2;

rea.mesh.nelg        = length(GLNew2);
rea.mesh.groupno     = rea.mesh.groupno(IndT);

% Sort them again
[GLNew2 I] = sort(GLNew2);
rea.mesh.xc  = rea.mesh.xc(:,I);
rea.mesh.yc  = rea.mesh.yc(:,I);
rea.mesh.cbc = rea.mesh.cbc(:,I);
rea.mesh.globalno    = GLNew2;
rea.mesh.groupno     = rea.mesh.groupno(I);

disp(['Total No. of Elements: ', num2str(rea.mesh.nelg)])

CheckConnectivity2D(rea.mesh)

CreateVTKMesh(rea.mesh) 

