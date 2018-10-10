function status = Nek_WriteRea(rea)

   fname = 'new.rea';
   wdsize = 'float32';   
%  Open file
   endian= 'le';
   [fid,message] = fopen(fname,'w+');
   
   if fid == -1
     disp(message) 
     return 
   end
 
   disp('Writing rea file...')

   nelg=rea.mesh.nelg;
   ndim=rea.mesh.ndim;
   ver=rea.nekver;

   WriteNekhdr(fid,ver,ndim);
   [pard, parc]=Nek_ParamDiscriptions;
   WriteParams(fid,rea,pard,parc); 

%  Skipping Passive scalar data. Needs to be coded in
   WriteNoPShdr(fid);

   nlogic=rea.nlogical;
   WriteLogicalSwitches(fid,rea,nlogic); 

   Nek_WriteReaMesh(rea.mesh,fid);   

   status = fclose(fid);


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
function WriteNoPShdr(fid)

%     This needs to be modified in a proper implementation

      hdr='Lines of passive scalar data follows 2 CONDUCT; 2RHOCP';
      nps = 0;
      space2=blanks(2);
      space5=blanks(5);
      fprintf(fid,'%s%2i%s%s\n',space5,nps,space2,hdr);

end   % function 
%----------------------------------------------------------------------
function WriteLogicalSwitches(fid,rea,nlogic)


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

%%----------------------------------------------------------------------

