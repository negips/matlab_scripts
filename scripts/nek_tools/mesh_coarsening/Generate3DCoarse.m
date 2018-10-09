function mesh3d = Generate3DCoarse(mesh2d,nz0,Lz,ifperiodic,ifvtk) 


      nlayers = length(mesh2d.layerindex);

      [mesh3d] = Generate3D(mesh2d,nlayers,nz0,Lz,ifperiodic);
      [~, nel]=size(mesh3d.xc);
      
      ndim=3;

      if (ifvtk)
        xvtk = zeros(2^ndim,nel);
        yvtk = zeros(2^ndim,nel);
        zvtk = zeros(2^ndim,nel);
        
        polydata = [];
        
        GLLind = [1 2 4 3 5 6 8 7];
        
        for i=1:nel
        
          xgll(:,i) = mesh3d.xc(GLLind,i);
          ygll(:,i) = mesh3d.yc(GLLind,i);
          zgll(:,i) = mesh3d.zc(GLLind,i);
        
          xvtk(:,i) = mesh3d.xc(:,i);
          yvtk(:,i) = mesh3d.yc(:,i);
          zvtk(:,i) = mesh3d.zc(:,i);
        
          p0 = (2^ndim)*(i-1);
        
          f1 = [0 1 5 4] + 1 + p0;
          f2 = [1 5 6 2] + 1 + p0;
          f3 = [2 6 7 3] + 1 + p0;
          f4 = [3 7 4 0] + 1 + p0;
          f5 = [0 1 2 3] + 1 + p0;
          f6 = [4 5 6 7] + 1 + p0;
        
          polydata = [polydata; f1; f2; f3; f4; f5; f6];
        
        end  
        
        vfname = 'mesh3d.vtk';
        vtkwrite(vfname,'polydata','hexahedron',xvtk,yvtk,zvtk,polydata,'binary')
      end   % ifvtk

end   % function
