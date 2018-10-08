function [mesh3d] = Build3DConnectivity(mesh3d,mesh2d)


   zlines = mesh2d.nelg;
   nz0 = length(mesh3d.LayerGEl{1});      % Maximum elelements in the span in layer1
   nfaces3 = 6;
   nfaces2 = 4;

   for il = 1:zlines         % 2d element number

     bcf1 = mesh2d.cbc(1,e2d).connectsto;    % bc on 2d face 1
     bcf2 = mesh2d.cbc(2,e2d).connectsto;    % bc on 2d face 2
     bcf3 = mesh2d.cbc(3,e2d).connectsto;    % bc on 2d face 3
     bcf4 = mesh2d.cbc(4,e2d).connectsto;    % bc on 2d face 4
        
     glnos = mesh3d.LayerGEl{e2d};
     nz=lnegth(glnos);
     cbc1 = [];

      
     for j=1:nz
       ieg = glnos(j);
%      Setting face boundary conditions
       for k=1:nfaces3
         cbc1(k).bc = mesh3d.EF{k,ieg};
       end  

%      Find connecting element
        
       
     end


end   % function

%----------------------------------------------------------------------

function ConnectivityS(mesh3d,mesh2d,e2d,j,k)



