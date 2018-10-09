function CheckConnectivity2D(mesh2d)

  nel    = length(mesh2d.globalno);
  disp(['Checking Element connectivity information for ' num2str(nel) ' elements'])

  nfaces = 4;
  ierr = 0;
  for i=1:nel
        
    elerr=0; 
    for j=1:nfaces
       bc = mesh2d.cbc(j,i).bc;
       if ~strcmpi(bc,'E  ')
%        Boundary condition
%        skip
         continue
       end
       ce = mesh2d.cbc(j,i).connectsto;
       of = mesh2d.cbc(j,i).onface;
       rev_el = mesh2d.cbc(of,ce).connectsto;
       if rev_el~=i
         disp(['Inconsistent Connectivity: ',num2str(j),' ', num2str(i), ' ', num2str(of), ' ', num2str(ce), ' ', num2str(rev_el)])
         ierr = ierr+1;
         elerr=elerr+1;
       end 
    end     %j=1:nfaces

    if elerr>0
      fig10=figure(10);
      Plot2DElement(mesh2d,i,fig10)
    end  

 end  % i=1:nel

 disp([num2str(ierr) ' inconsistencies found'])

end   % function
%---------------------------------------------------------------------- 

