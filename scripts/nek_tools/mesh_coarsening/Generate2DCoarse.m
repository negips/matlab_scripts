function [rea2d] = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk)
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

%       If first few layers are smaller than the rest        
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
%         ifc(j) = CoarsenSaab(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenSaab600k(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenSaab2d(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenNaca2d(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenNaca2d_2(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenNaca77k1(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenSkern(LX,LY,j,i,iflocked);
%         ifc(j) = CoarsenNaca77k_5_10_fine2(LX,LY,j,i,iflocked);
         ifc(j) = CoarsenNacatrunc(LX,LY,j,i,iflocked);

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

      rea2d.casename          = rea.casename;
      rea2d.nekver            = rea.nekver;
      rea2d.nparams           = rea.nparams;
      rea2d.Nlogical          = rea.Nlogical;
      rea2d.npscal            = rea.npscal;
      rea2d.ifflow            = rea.ifflow;
      rea2d.ifheat            = rea.ifheat;
      rea2d.param             = rea.param;
      rea2d.logical           = rea.logical;
      rea2d.mesh              = mesh2d;
      rea2d.Nrestart          = rea.Nrestart;
      rea2d.rstFiles          = rea.rstFiles;
      rea2d.rstOptions        = rea.rstOptions;
      rea2d.Nic               = rea.Nic;
      rea2d.initialconditions = rea.initialconditions;
      rea2d.Ndriveforce       = rea.Ndriveforce;
      rea2d.driveforce        = rea.driveforce;
      rea2d.Nvplines          = rea.Nvplines;
      rea2d.Npackets          = rea.Npackets;
      rea2d.datapacket        = rea.datapacket;
      rea2d.Nhist             = rea.Nhist;
      rea2d.history           = rea.history;
      rea2d.Noutspec          = rea.Noutspec;
      rea2d.outputspec        = rea.outputspec;
      rea2d.Nobjects          = rea.Nobjects;
      rea2d.objects           = rea.objects;
      
      if ifvtk
        CreateVTKMesh(rea2d.mesh) 
      end  

end   % function


 



