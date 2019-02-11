function CreateVTKMesh(mesh)

      ndim=mesh.ndim;
      if ndim==2
        vfname='mesh2d.vtk';
        polytype='tetrahedron';
      elseif ndim==3
        vfname='mesh3d.vtk';
        polytype='hexahedron';
      else
        disp(['Not implemented for ndim=', nu2str(ndim)])
        return
      end

      disp(['Writing ' vfname])

      nel=mesh.nelg;
      xvtk = zeros(2^ndim,nel);
      yvtk = zeros(2^ndim,nel);
      zvtk = zeros(2^ndim,nel);
      if ndim==2
        npoly=1;
      elseif ndim==3
        npoly=6;
      end  
      polydata = zeros(npoly*nel,4);

      ip = 0;
      for i=1:nel

        xvtk(:,i) = mesh.xc(:,i);
        yvtk(:,i) = mesh.yc(:,i);
        if ndim==3
          zvtk(:,i) = mesh.zc(:,i);
        end
      
        p0 = (2^ndim)*(i-1);
        if ndim==3      
          f1 = [0 1 5 4] + 1 + p0;
          f2 = [1 5 6 2] + 1 + p0;
          f3 = [2 6 7 3] + 1 + p0;
          f4 = [3 7 4 0] + 1 + p0;
          f5 = [0 1 2 3] + 1 + p0;
          f6 = [4 5 6 7] + 1 + p0;
%          polydata = [polydata; f1; f2; f3; f4; f5; f6];
          polydata(ip+1:ip+npoly,:) = [f1; f2; f3; f4; f5; f6];
          ip = ip+npoly;
        else
          f1 = [0 1 2 3] + 1 + p0;
%          polydata = [polydata; f1];
          ip=ip+1;
          polydata(ip,:) = f1;
        end  
      end
      
      xvtk=xvtk(:);
      yvtk=yvtk(:);
      zvtk=zvtk(:);
      
      vtkwrite(vfname,'polydata',polytype,xvtk,yvtk,zvtk,polydata,'binary')

end   % function
