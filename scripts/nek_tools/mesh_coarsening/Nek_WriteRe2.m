function status = Nek_WriteRe2(mesh)

   fname = 'mesh.re2';
   wdsize = 'float32';   
%  Open file
   endian= 'le';
   [fid,message] = fopen(fname,'w+',['ieee-' endian]);
   
   if fid == -1
     disp(message) 
     return 
   end
 
   disp('Writing re2 file...')

   nelg=mesh.nelg;
   ndim=mesh.ndim;

   hdr = sprintf('#v002%9i%3i%9i hdr',nelg,ndim,nelg);
   l=length(hdr);
   for i=l+1:80
      hdr = [hdr ' '];
   end
   length(hdr)
   s = uint8(hdr);   
   fwrite(fid,s,'char') 

%  Endian tag
   etag = 6.54321;
   fwrite(fid,etag,wdsize);

   fclose(fid)







end   % function

