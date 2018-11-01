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

      if ifvtk
        CreateVTKMesh(rea3d.mesh) 
      end  
     

end   % function
