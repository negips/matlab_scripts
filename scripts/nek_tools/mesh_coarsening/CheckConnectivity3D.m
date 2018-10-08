function CheckConnectivity3D(mesh3d)


  nel    = length(mesh3d.globalno);
  disp(['Checking Element connectivity information for ' num2str(nel) ' elements'])

  nfaces = 6;
  ierr = 0;
  for i=1:nel
        
    elerr=0; 
    for j=1:nfaces
       bc = mesh3d.cbc(j,i).bc;
       if ~strcmpi(bc,'E  ') && ~strcmpi(bc,'P  ')
%        Boundary condition
%        skip
         continue
       end
       ce = mesh3d.cbc(j,i).connectsto;
       of = mesh3d.cbc(j,i).onface;
       rev_el = mesh3d.cbc(of,ce).connectsto;
       if rev_el~=i
         disp(['Inconsistent Connectivity ',num2str(j),' ', num2str(i), ' ', num2str(of), ' ', num2str(ce), ' ', num2str(rev_el)])
         ierr = ierr+1;
         elerr=elerr+1;
       end 
    end     %j=1:nfaces

    if elerr>0
      fig10=figure(10);
      Plot3DElement(mesh3d,i,fig10)
    end  

 end  % i=1:nel

 disp([num2str(ierr) ' inconsistencies found'])

end   % function
%---------------------------------------------------------------------- 

