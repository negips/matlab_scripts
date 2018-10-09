function [mesh2d] = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk)
%     First attempts at coarsening the mesh

%           Element arrangement
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
      
      ifplot = 0;
      
      disp(['Total Number of Elements: ', num2str(rea.mesh.nelg)])
      
      nlayers = length(LayerE);
      
      for i=1:nlayers
        j=nlayers-i+1;
        NewE{i}=LayerE{j};
        NewX{i}=LayerX{j};
        NewY{i}=LayerY{j};
        NewBC{i}=LayerBC{j};
        NewCEl{i}=LayerCEl{j};
        NewMeshC(i)=MeshC(j);
      
%       Define an element type 
        l1=length(NewE{i});
        for j=1:l1
          NewET{i}{j}='s';
        end  
      
      end
      
      NewCoF = ConnectedOnFace(NewE,NewCEl,NewMeshC);
      
%     Coarsen Layer by Layer
%     Define aspect ratio as 'O' face lengths to 'V' face lengths
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
      
          continue
        end  
      
        l1 =length(LX);
        if i==layer_start
          iflocked=zeros(l1,1);
        end  
        cmap = jet(l1);
        l2 = length(NewX{i});
        ifc = zeros(l2,1); 
        for j=1:l1
            
%         ifc(j) = CoarsenCriteria(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenKDJ(LX,LY,j,i,iflocked);
         ifc(j) = CoarsenSaab(LX,LY,j,i,iflocked);
      
        end

%       Just ensure we don't try to coarsen last layer      
        if i==nlayers
          ifc(:)=0;
        end  
      
        c_ind = find(ifc);
        nc = length(c_ind);
        if nc>0
          disp(['Coarsen ' num2str(nc) ' Elements in layer ' num2str(i)])
        end  
      
        ifplot =0;
%       Coarsen layer
        [LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,NewET]=Coarsen2DLayer(LE,LX,LY,ifc,NewX,NewY,NewCEl,NewCoF,NewMeshC,NewET,i,fig2,ifplot);
%       Modify all subsequent Layers
        [NewE, NewX, NewY, NewBC, NewCEl, NewCoF, NewET, iflocked]=Create2DLayers(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,NewET,i,ifc,fig2,ifplot);
        nels_layer = length(NewX{i});
        new_nelg = new_nelg+nels_layer;
      end
      disp(['New Number of Elements: ', num2str(new_nelg)])
      
%     Not done for multiple curvature definitions right now
      mesh2d = ReOrderElements(NewE,NewX,NewY,NewBC,NewCEl,NewCoF,NewET,rea.mesh,curvedef); 
      CheckConnectivity2D(mesh2d)
      
      if ifvtk      
        polydata = [];
        nel=mesh2d.nelg;
        for i=1:nel
          xvtk(:,i) = mesh2d.xc(:,i);
          yvtk(:,i) = mesh2d.yc(:,i);
        
          p0 = (2^2)*(i-1);
          f1 = [0 1 2 3] + 1 + p0;
          polydata = [polydata; f1];
        end  
        
        xvtk=xvtk(:);
        yvtk=yvtk(:);
        zvtk=0*xvtk;
        
        vfname = 'mesh2d.vtk';
        vtkwrite(vfname,'polydata','tetrahedron',xvtk,yvtk,zvtk,polydata,'binary')
      end  

end   % function


 



