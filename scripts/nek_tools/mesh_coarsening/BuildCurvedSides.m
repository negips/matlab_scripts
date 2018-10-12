function mesh3d = BuildCurvedSides(mesh3d,mesh2d,Lz)

%  Also add curved element information
%  I assume the curved elements have not been modified.
%  Otherwise the code needs a lot of modifications
   
   ncurves = mesh2d.Ncurve;
   curveparams = [];
   curveedge   = [];
   curveieg    = [];
   ctype       = [];   
   nc3d  = 0;
   nedge = 0;
   for i=1:ncurves
     el2d = mesh2d.curveieg(i);
     f2d  = mesh2d.curveedge(i);
     ct   = mesh2d.curvetype{i};
     if ~strcmpi(ct,'m')
       disp(['Not implemented for Curve Type ', ct, '. Skipping curved element'])
       continue
     end

     nel3d=length(mesh3d.LayerGEl{el2d});
     dz = Lz/nel3d;
     for j=1:nel3d
       nc3d = nc3d+1;
       for k=0:1
         nedge   = nedge+1;
         cp = zeros(5,1);
         cieg    =  mesh3d.LayerGEl{el2d}(j);     % 3d Element number
         iedge   =  mesh2d.curveedge(i) + k*4;    % edge number
         cp(1)   =  mesh2d.curveparams(1,i);
         cp(2)   =  mesh2d.curveparams(2,i);
         cp(3)   =  ((j-1)*2 + k)*dz;
         cp(4)   = 0.;
         cp(5)   = 0.;

         curveieg  = [curveieg  cieg];
         curveedge = [curveedge iedge];
         curveparams = [curveparams cp];
         ctype{nedge} = 'm';
       end  % k=0:1
     end    % j=1:nel3d
   end      % i=1:ncurves 
      
   mesh3d.Ncurve      = nedge;
   mesh3d.curveedge   = curveedge;
   mesh3d.curveieg    = curveieg;
   mesh3d.curveparams = curveparams;
   mesh3d.curvetype   = ctype;

end   % function
%---------------------------------------------------------------------- 
