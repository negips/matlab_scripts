function rea3d = Generate3DCoarse(rea2d,nz0,Lz,ifperiodic,ifvtk) 

      nlayers = length(rea2d.mesh.layerindex);
      [mesh3d] = Generate3D(rea2d.mesh,nlayers,nz0,Lz,ifperiodic);
      nel=mesh3d.nelg;

      rea3d.casename          = rea2d.casename;
      rea3d.nekver            = rea2d.nekver;
      rea3d.nparams           = rea2d.nparams;
      rea3d.Nlogical          = rea2d.Nlogical;
      rea3d.npscal            = rea2d.npscal;
      rea3d.ifflow            = rea2d.ifflow;
      rea3d.ifheat            = rea2d.ifheat;
      rea3d.param             = rea2d.param;
      rea3d.logical           = rea2d.logical;
      rea3d.mesh              = mesh3d;
      rea3d.Nrestart          = rea2d.Nrestart;
      rea3d.rstFiles          = rea2d.rstFiles;
      rea3d.rstOptions        = rea2d.rstOptions;
      rea3d.Nic               = rea2d.Nic;
      rea3d.initialconditions = rea2d.initialconditions;
      rea3d.Ndriveforce       = rea2d.Ndriveforce;
      rea3d.driveforce        = rea2d.driveforce;
      rea3d.Nvplines          = rea2d.Nvplines;
      rea3d.Npackets          = rea2d.Npackets;
      rea3d.datapacket        = rea2d.datapacket;
      rea3d.Nhist             = rea2d.Nhist;
      rea3d.history           = rea2d.history;
      rea3d.Noutspec          = rea2d.Noutspec;
      rea3d.outputspec        = rea2d.outputspec;
      rea3d.Nobjects          = rea2d.Nobjects;
      rea3d.objects           = rea2d.objects;

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
