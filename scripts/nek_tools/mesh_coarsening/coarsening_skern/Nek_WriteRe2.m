function Nek_WriteRe2(mesh)

      disp('Writing re2')

      endian= 'le';
%     Double precision   
      rt = 'float64';

      nelg=mesh.nelg;
      ndim=mesh.ndim;
      nelv=nelg;

      if ndim==3
        fname2='new3d.re2';
      else
        fname2='new2d.re2';
      end  

      [fre2,message] = fopen(fname2,'w+',['ieee-' endian]);
      if fre2 == -1
        disp(message) 
        return 
      end


%     header   
      hdr = sprintf('#v002%9i%3i%9i%s',nelg,ndim,nelv,' hdr');

      l1=length(hdr);
      hdrlen=80;
      if l1<hdrlen
        space=blanks(hdrlen-l1);
        hdr=[hdr space];
      else
        disp(['Truncating header to ' num2str(hdrlen) ' characters']);
        hdr=hdr(1:hdrlen);
      end
     
      fwrite(fre2,hdr,'char');
%     endian discriminator   
      etag = 6.54321;
      fwrite(fre2,etag,'float32');

%     Mesh data   
      igroup = 0;
      for e=1:nelg 
        WriteXYZ(fre2,e,igroup,mesh,ndim,rt);
      end

%     Curved sides   
      ncurve=mesh.Ncurve;
      fwrite(fre2,ncurve,rt);       % No of curves 
      for e=1:ncurve      
        WriteCurveData(fre2,mesh,e,rt)
      end   

%     Fluid BCs
      nbcs = FindBCs(mesh,ndim,nelg);           % no of BCs
      fwrite(fre2,nbcs,rt); 
      for e=1:nelg
         WriteFluidBC(fre2,mesh,e,ndim,rt) 
      end

      fclose(fre2);

end   % function
%---------------------------------------------------------------------- 
function WriteXYZ(fid,e,igroup,mesh,ndim,rt)

      fwrite(fid,double(igroup),rt);
      if ndim==2
        xc=mesh.xc(:,e)';
        fwrite(fid,xc,rt);
        yc=mesh.yc(:,e)';
        fwrite(fid,yc,rt);
      else
        xc=mesh.xc(:,e);
        fwrite(fid,xc,rt);
        yc=mesh.yc(:,e);
        fwrite(fid,yc,rt);
        zc=mesh.zc(:,e);
        fwrite(fid,zc,rt);
      end   

end   % function
%----------------------------------------------------------------------
function  WriteCurveData(fid,mesh,e,rt)

      ieg=double(mesh.curveieg(e));
      edge=double(mesh.curveedge(e));
      cp1=mesh.curveparams(1,e);
      cp2=mesh.curveparams(2,e);
      cp3=mesh.curveparams(3,e);
      cp4=mesh.curveparams(4,e);
      cp5=mesh.curveparams(5,e);
      buf = [ieg edge cp1 cp2 cp3 cp4 cp5];
      fwrite(fid,buf,rt);

      ct = mesh.curvetype{e};
      l1=length(ct);

      ct1=ct(1);              % Buffer length in Nek is 1 real
      fwrite(fid,ct1,rt);

end   % function
%---------------------------------------------------------------------- 
function nbcs = FindBCs(mesh,ndim,nelg)
      
      nbcs = 0;
      nfaces=2*ndim;
      for i=1:nelg
        for j=1:nfaces
          bc=mesh.cbc(j,i).bc;
          if ~strcmpi(bc,'E  ')
             nbcs=nbcs+1;
          end
        end
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteFluidBC(fid,mesh,e,ndim,rt)

      nfaces=2*ndim;
      for j=1:nfaces
        bc=mesh.cbc(j,e).bc;

        if strcmpi(bc,'E  ')
          continue
        end  

        c2=double(mesh.cbc(j,e).connectsto);
        of=double(mesh.cbc(j,e).onface);
        if isfield(mesh.cbc(j,e),'param3');
          p3=mesh.cbc(j,e).param3;
        else
          p3 = 0.;
        end  
        if isfield(mesh.cbc(j,e),'param4');
          p4=mesh.cbc(j,e).param4;
        else
          p4 = 0.;
        end  
        if isfield(mesh.cbc(j,e),'param5');
          p5=mesh.cbc(j,e).param5;
        else
          p5 = 0.;
        end 

        buf = [double(e) double(j) c2 of p3 p4 p5];
        fwrite(fid,buf,rt);

        l1=length(bc);
        CL=8;           % Buffer length in nek for BCs
        if l1>CL
          bc=bc(1:CL);
        elseif l1<CL
          s=blanks(CL-l1);
          bc = [bc s];
        end 

        fwrite(fid,bc,'char*1');
      end

end   % function
%----------------------------------------------------------------------

