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
%%----------------------------------------------------------------------
function WriteParams(fid,rea,pard,parc)

      hdr=' PARAMETERS FOLLOW';
      space5=blanks(5);             % arbitrary
      nparams=rea.nparams;
      fprintf(fid,'%s%5i%s\n',space5,nparams,hdr);

      fmt='%s%14.6e%s\n';
      disfmt=repmat('%s',1,6);

      space2=blanks(2);
      space4=blanks(4);
      for i=1:85 %nparams
        pn = sprintf('P%3.3i',i);
        l1=length(pard{i});
        spacen=blanks(16-l1);
        description=sprintf(disfmt,space2,pn,space2,pard{i},spacen,parc{i});

        par=rea.param(i);
        fprintf(fid,fmt,space2,par,description);
      end  


end   % function
%---------------------------------------------------------------------- 

%---------------------------------------------------------------------- 

%function WriteMeshhdr(fid,ndim,nelg)
%
%      if ndim==3
%        hdr=' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8';
%        fprintf(fid,'%s\n',hdr);
%      else
%        hdr=' **MESH DATA** 2 lines are X,Y. Columns corners 1-4';
%        fprintf(fid,'%s\n',hdr);
%      end  
%
%      space11=blanks(11);
%      hdr='NEL,NDIM,NELV';
%      fprintf(fid,'%12i%3i%12i%s%s\n',nelg,ndim,nelg,space11,hdr); 
%
%end   % function
%%----------------------------------------------------------------------

