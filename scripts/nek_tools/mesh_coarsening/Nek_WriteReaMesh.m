function Nek_WriteReaMesh(mesh,fid)

   disp('Writing mesh')

   nelg=mesh.nelg;
   ndim=mesh.ndim;

%  Xfac,Yfac,Xzero,YZero
   WritePrenekhdr(fid,mesh)

%  Mesh hdr
   WriteMeshhdr(fid,ndim,nelg);

%  Mesh data   
   if ndim==2
     zlev = 1;
     a = ' ';
     igroup = 0;
     for e=1:nelg 
       WriteElementhdr(fid,e,zlev,a,igroup);
       WriteXYZRea(fid,e,mesh,ndim);
     end
   else % 3d

   end       

%  Curved sides   
   ncurve=mesh.ncurve;
   ncurve=0.
   WriteCurvehdr(fid,ncurve)   
   for e=1:ncurve      
     WriteCurveData(fid,mesh,e,nelg)
   end   

%  Fluid BCs
   WriteFluidBChdr(fid)
   for e=1:nelg
      WriteFluidBC(fid,mesh,e,ndim,nelg)
   end

%  Have not implemented thermal boundary conditions
   WriteNoThermalhdr(fid)             

end   % function
%---------------------------------------------------------------------- 
function WritePrenekhdr(fid,mesh)

      space5 = blanks(5);
      space4 = blanks(4);

      hdr='XFAC,YFAC,XZERO,YZERO';
         
      fprintf(fid,'%10f%s%10f%s%10f%s%10f%s%s\n',mesh.xfac,space4,mesh.yfac,space4,mesh.xzero,space4,mesh.yzero,space5,hdr);

end   % function 
%----------------------------------------------------------------------  

function WriteMeshhdr(fid,ndim,nelg)

      if ndim==3
        hdr=' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8';
        fprintf(fid,'%s\n',hdr);
      else
        hdr=' **MESH DATA** 2 lines are X,Y. Columns corners 1-4';
        fprintf(fid,'%s\n',hdr);
      end  

      space11=blanks(11);
      hdr='NEL,NDIM,NELV';
      fprintf(fid,'%12i%3i%12i%s%s\n',nelg,ndim,nelg,space11,hdr); 

end   % function
%----------------------------------------------------------------------
function WriteElementhdr(fid,e,zlev,a,igroup)

      space5 = blanks(5);
      space4 = blanks(4);
         
      fprintf(fid,'%sELEMENT%12i [%5i%1s]%sGROUP%5i\n',space5,e,zlev,a,space4,igroup);

end   % function 
%----------------------------------------------------------------------  
function WriteXYZRea(fid,e,mesh,ndim)

      fmt = repmat('%14.6e',1,4);
      fmt = [fmt '\n'];    
         
      if ndim==2   
        fprintf(fid,fmt,mesh.xc(:,e)); 
        fprintf(fid,fmt,mesh.yc(:,e));
      else
        ind=1:4; 
        fprintf(fid,fmt,mesh.xc(ind,e)); 
        fprintf(fid,fmt,mesh.yc(ind,e));
        fprintf(fid,fmt,mesh.zc(ind,e));
        ind=5:8; 
        fprintf(fid,fmt,mesh.xc(ind,e)); 
        fprintf(fid,fmt,mesh.yc(ind,e));
        fprintf(fid,fmt,mesh.zc(ind,e));
      end   

end   % function
%----------------------------------------------------------------------
function WriteCurvehdr(fid,ncurve)

      hdr = '  ***** CURVED SIDE DATA *****';
      fprintf(fid,'%s\n',hdr);
      hdr = 'Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE';
      space2=blanks(2);
      fprintf(fid,'%12i%s%s\n',ncurve,space2,hdr);

end   % function 
%----------------------------------------------------------------------
function  WriteCurveData(fid,mesh,e,nelg)

      datafmt = repmat('%14.6e',1,5);
      space1=blanks(1);
      if nelg<1000
        fmt=['%3i%3i',datafmt,space1,'%1s\n'];
      elseif nelg<10^6
        fmt=['%2i%6i',datafmt,space1,'%1s\n'];
      else
        fmt=['%2i%12i',datafmt,space1,'%1s\n'];
      end

      ieg=mesh.curveieg(e);
      edge=mesh.curveface(e);
      cp1=mesh.curveparams(1);
      cp2=mesh.curveparams(2);
      cp3=mesh.curveparams(3);
      cp4=mesh.curveparams(4);
      cp5=mesh.curveparams(5);
      ct =mesh.curvetype{e};

      fprintf(fid,fmt,edge,ieg,cp1,cp2,cp3,cp4,cp5,ct);

%     60 format(i3,i3,5g14.6,1x,a1)
%     61 format(i2,i6,5g14.6,1x,a1)
%     62 format(i2,i12,5g14.6,1x,a1)

end   % function
%---------------------------------------------------------------------- 
function WriteFluidBChdr(fid)

      hdr = '  ***** BOUNDARY CONDITIONS ****';
      fprintf(fid,'%s\n',hdr);
      hdr = '  ***** FLUID   BOUNDARY CONDITIONS *****';
      fprintf(fid,'%s\n',hdr);

end   % function 
%---------------------------------------------------------------------- 
function WriteFluidBC(fid,mesh,e,ndim,nelg)

%     Either there is a bug for large element numbers or 
%     there is an inherent ordering of elements for BCs.

      datafmt = repmat('%14.6e',1,5);
      space1  = blanks(1);

      fmt1=0;

      if nelg<1000      
        fmt = [space1,'%3s%3i%3i',datafmt,'\n'];
        fmt1=1;
      elseif nelg<10^5
        fmt = [space1,'%3s%6i%1i',datafmt,'\n'];
        fmt1=1;
      elseif nelg<10^6
        fmt = [space1,'%3s%5i%1i',datafmt,'\n'];
        fmt1=0;
      else
        datafmt = repmat('%18.11',1,5);
        fmt = [space1,'%3s%12i',datafmt,'\n'];
        fmt1=0;
      end

      nfaces=2*ndim;
      for j=1:nfaces
        bc=mesh.cbc(j,e).bc;
        c2=mesh.cbc(j,e).connectsto;
        of=mesh.cbc(j,e).onface;
        p1=mesh.cbc(j,e).param1;
        p2=mesh.cbc(j,e).param2;
        p3=mesh.cbc(j,e).param3;

        if fmt1
          fprintf(fid,fmt,bc,e,j,c2,of,p1,p2,p3);
        else
          fprintf(fid,fmt,bc,e,c2,of,p1,p2,p3);
        end
      end

%     20 FORMAT(1x,A3,2I3,5G14.6)
%     21 FORMAT(1x,A3,i5,i1,5G14.6)
%     22 FORMAT(1x,A3,i6,5G14.6)
%     23 FORMAT(1x,A3,i12,5G18.11)

end   % function
%----------------------------------------------------------------------
function WriteNoThermalhdr(fid)

      hdr = '  ***** NO THERMAL BOUNDARY CONDITIONS *****';
      fprintf(fid,'%s\n',hdr);

end   % function 
%---------------------------------------------------------------------- 

