function status = Nek_WriteRea(rea,ifre2)

   nelg=rea.mesh.nelg;
   nelv=nelg;
   ndim=rea.mesh.ndim;
   ver=rea.nekver;

   if ndim==3
     fname = 'new3d.rea';
   else  
     fname = 'new2d.rea';
   end  

   disp(['Writing rea file ' fname])

   wdsize = 'float64';   
%  Open file
   [frea,message] = fopen(fname,'w+');
   
   if frea == -1
     disp(message) 
     return 
   end



   WriteNekhdr(frea,ver,ndim);
   [pard, parc]=Nek_ParamDiscriptions;
   WriteParams(frea,rea,pard,parc); 

%  Skipping Passive scalar data. Needs to be coded in
   WriteNoPShdr(frea);

%  Logical Switches   
   WriteLogicalSwitches(frea,rea); 

%  Xfac,Yfac,Xzero,YZero
   WritePrenekhdr(frea,rea.mesh)

%  Mesh hdr
   WriteMeshhdr(frea,ndim,nelg,nelv,ifre2);

%  Mesh data
   if ifre2
     Nek_WriteRe2(rea.mesh);
   else            
     Nek_WriteReaMesh(frea,rea.mesh);   
   end

%  Restart
   WriteRestart(frea,rea)

%  Initial
   WriteInitialCond(frea,rea)

%  Drive Force
   WriteDriveForce(frea,rea)

%  Drive Force
   WriteVarProp(frea,rea)

%  History and Integral
   WriteHist(frea,rea)

%  Output specification
   WriteOutputSpec(frea,rea)

%  Object specification
   WriteObjectSpec(frea,rea)
  
   status = fclose(frea);


end   % function
%---------------------------------------------------------------------- 
function WriteNekhdr(fid,ver,ndim)

      hdr='****** PARAMETERS ******';
      fprintf(fid,'%s\n',hdr);
      
      space5=blanks(5);
      hdr='NEKTON VERSION';
      fprintf(fid,'%12.6f%s%s\n',ver,space5,hdr);

      space2=blanks(2);
      hdr=' DIMENSIONAL RUN';
      fprintf(fid,'%s%2i%s\n',space2,ndim,hdr);

end   % function 
%----------------------------------------------------------------------
function WriteParams(fid,rea,pard,parc)

      hdr=' PARAMETERS FOLLOW';
      space5=blanks(5);             % arbitrary
      nparams=rea.nparams;
      fprintf(fid,'%s%5i%s\n',space5,nparams,hdr);

      fmt='%s%14.6e%s\n';
      disfmt=repmat('%s',1,6);

      space2=blanks(2);
      space4=blanks(4);
      for i=1:nparams
        pn = sprintf('P%3.3i',i);
        l1=length(pard{i});
        spacen=blanks(16-l1);
        description=sprintf(disfmt,space2,pn,space2,pard{i},spacen,parc{i});

        par=rea.param(i);
        fprintf(fid,fmt,space2,par,description);
      end  


end   % function
%---------------------------------------------------------------------- 
function WritePrenekhdr(fid,mesh)

      space5 = blanks(5);
      space4 = blanks(4);
      hdr='XFAC,YFAC,XZERO,YZERO';

      datafmt = repmat('%10f%s',1,4);
      fmt     = [datafmt '%s\n'];
      fprintf(fid,fmt,mesh.xfac,space4,mesh.yfac,space4,mesh.xzero,space4,mesh.yzero,space5,hdr);

end   % function
%---------------------------------------------------------------------- 
function WriteMeshhdr(fid,ndim,nelg,nelv,ifre2)

      if ndim==3
        hdr=' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8';
        fprintf(fid,'%s\n',hdr);
      else
        hdr=' **MESH DATA** 2 lines are X,Y. Columns corners 1-4';
        fprintf(fid,'%s\n',hdr);
      end  

      space11=blanks(11);
      hdr='NEL,NDIM,NELV';
      if ifre2
        nelg=-nelg;
      end  
      fprintf(fid,'%12i%3i%12i%s%s\n',nelg,ndim,nelv,space11,hdr); 

end   % function
%----------------------------------------------------------------------
function WriteNoPShdr(fid)

%     This needs to be modified in a proper implementation

      hdr='Lines of passive scalar data follows 2 CONDUCT; 2RHOCP';
      nps = 0;
      space2=blanks(2);
      space5=blanks(5);
      fprintf(fid,'%s%2i%s%s\n',space5,nps,space2,hdr);

end   % function 
%----------------------------------------------------------------------
function WriteLogicalSwitches(fid,rea)

      nlogic=rea.Nlogical;

%     header
      space2=blanks(2);
      space5=blanks(5);
      hdr='LOGICAL SWITCHES FOLLOW';
      fprintf(fid,'%s%3i%s%s\n',space5,nlogic,space2,hdr);

%     Switches      
      for i=1:nlogic
        flag=rea.logical{i,2};
        val =rea.logical{i,1};            
        l1=length(val);
        if l1==1
          fprintf(fid,'%s%s%s%s\n',space2,val,space5,flag);
        else
          val2=[];
          for j=1:l1
            val2 = [val2 val(j) ' '];
          end
          val2(end)=[]; % remove last space
          fprintf(fid,'%s%s%s%s\n',space2,val2,space2,flag);
        end  
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteRestart(fid,rea)

      nrestart=rea.Nrestart;

      space5 = blanks(5);
      space1 = blanks(1);
      hdr='PRESOLVE/RESTART OPTIONS';

      fmt     = '%s%3i%s%s\n';
      fprintf(fid,fmt,space5,nrestart,space1,hdr);

      fmt     = '%s%s%s\n';
      for i=1:nrestart
        fprintf(fid,fmt,rea.rstFiles{i},space1,rea.rstOptions{i});
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteInitialCond(fid,rea)

      nic=rea.Nic;

      space5 = blanks(5);
      space1 = blanks(1);
      hdr='INITIAL CONDITIONS';

      fmt     = '%s%3i%s%s\n';
      fprintf(fid,fmt,space5,nic,space1,hdr);

      fmt     = '%s\n';
      for i=1:nic
        fprintf(fid,fmt,rea.initialconditions{i});
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteDriveForce(fid,rea)

      ndrive=rea.Ndriveforce;

      hdr=' ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q ';
      fmt     = '%s\n';
      fprintf(fid,fmt,hdr);

      space5 = blanks(5);
      space1 = blanks(1);
      hdr='Lines of Drive force data follow';

      fmt     = '%s%3i%s%s\n';
      fprintf(fid,fmt,space5,ndrive,space1,hdr);

      fmt     = '%s\n';
      for i=1:ndrive
        fprintf(fid,fmt,rea.driveforce{i});
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteVarProp(fid,rea)

      hdr=' ***** VARIABLE PROPERTY DATA ***** Overrides Parameter data ';
      fmt     = '%s\n';
      fprintf(fid,fmt,hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      hdr='Lines follow';
      fmt     = '%s%3i%s%s\n';
      fprintf(fid,fmt,space5,rea.Nvplines,space1,hdr);

      hdr='PACKETS OF DATA FOLLOW';
      fmt     = '%s%3i%s%s\n';
      npackets = rea.Npackets;
      fprintf(fid,fmt,space5,npackets,space1,hdr);
     
      fmt     = '%s\n';
      for i=1:npackets
        fprintf(fid,fmt,rea.datapacket{i});
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteHist(fid,rea)

      nhist=rea.Nhist;

      hdr=' ***** HISTORY AND INTEGRAL DATA *****';
      fmt     = '%s\n';
      fprintf(fid,fmt,hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      hdr='POINTS. Hcode, I,J,H,IEL';
      fmt     = '%s%3i%s%s\n';
      fprintf(fid,fmt,space5,nhist,space1,hdr);
     
      fmt     = '%s\n';
      for i=1:nhist
        fprintf(fid,fmt,rea.history{i});
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteOutputSpec(fid,rea)

      noutspec=rea.Noutspec;

      hdr=' ***** OUTPUT FIELD SPECIFICATION *****';
      fmt     = '%s\n';
      fprintf(fid,fmt,hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      hdr='SPECIFICATIONS FOLLOW';
      fmt     = '%s%3i%s%s\n';
      fprintf(fid,fmt,space5,noutspec,space1,hdr);

      space2 = blanks(2);
      fmt     = '%s%s%s%s\n';
      for i=1:noutspec
        fprintf(fid,fmt,space2,rea.outputspec{i,1},space5,rea.outputspec{i,2});
      end  

end   % function
%---------------------------------------------------------------------- 
function WriteObjectSpec(fid,rea)

      nobject=rea.Nobjects;

      hdr=' ***** OBJECT SPECIFICATION *****';
      fmt     = '%s\n';
      fprintf(fid,fmt,hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      space2 = blanks(2);
      fmt     = '%s%3i%s%s\n';
      for i=1:nobject
        fprintf(fid,fmt,space2,rea.objects{i,1},space5,rea.objects{i,2});
      end  

end   % function
%---------------------------------------------------------------------- 

%    From genxyz.f in genbox:
%C   Preprocessor Corner notation:      Symmetric Corner notation:
%C
%C           4+-----+3    ^ s                    3+-----+4    ^ s
%C           /     /|     |                      /     /|     |
%C          /     / |     |                     /     / |     |
%C        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
%C         |     | /     /                     |     | /     /
%C         |     |/     /                      |     |/     /
%C        5+-----+6    t                      5+-----+6    t
%C
