function [mesh3d] = Build3DConnectivity(mesh3d,mesh2d)


   zlines = mesh2d.nelg;
   nz0 = length(mesh3d.LayerGEl{1});      % Maximum elelements in the span in layer1
   nfaces3 = 6;
   nfaces2 = 4;

   nlayers=length(mesh2d.layerindex);
   rsum_z = [];
   rsum = 0;
   e=0;
   for il=1:nlayers
     lind=mesh2d.layerindex{il};
     l_nel=length(lind);
     for j=1:l_nel
       e=e+1;
       nzl = length(mesh3d.LayerGEl{e});          % No of elements in this 'z' line
       rsum = rsum+nzl;
       rsum_z(e)=rsum;
     end
   end

%  Running sum of global element numbers 
   mesh3d.rsum_z = rsum_z; 

   cbc1 = [];
   e3d=0;
   nlayers=length(mesh2d.layerindex);
   for il = 1:nlayers
     disp(['Building Connectivity for layer ', num2str(il)])
     lind_2d=mesh2d.layerindex{il};
     for i = lind_2d         % 2d element number

       bcf1 = mesh2d.cbc(1,i).connectsto;    % bc on 2d face 1
       bcf2 = mesh2d.cbc(2,i).connectsto;    % bc on 2d face 2
       bcf3 = mesh2d.cbc(3,i).connectsto;    % bc on 2d face 3
       bcf4 = mesh2d.cbc(4,i).connectsto;    % bc on 2d face 4
          
       glnos = mesh3d.LayerGEl{i};
       nz=length(glnos);
        
       if mesh3d.ifzc(il)==0           % No spanwise coarsening

         for j=1:nz
           e3d=e3d+1;
           ieg = glnos(j);
%          Setting face boundary conditions
           for k=1:nfaces3
             cbc1(k).bc = mesh3d.EF{ieg}(k,:);
           end  
%          Find connecting element on face 5
           if strcmpi(cbc1(5).bc,'P  ')          % if its a periodic face
             cbc1(5).connectsto = ieg+nz-1;
             cbc1(5).onface     = 6;
           elseif ~strcmpi(cbc1(5).bc,'E  ')     % If its a different boundary condition
             cbc1(5).connectsto = 0;
             cbc1(5).onface     = 0;
           else
             cbc1(5).connectsto = ieg-1;
             cbc1(5).onface     = 6;
           end  
%          Find connecting element on face 6
           if strcmpi(cbc1(6).bc,'P  ')
             cbc1(6).connectsto = ieg-nz+1;
             cbc1(6).onface     = 5;
           elseif ~strcmpi(cbc1(6).bc,'E  ')
             cbc1(6).connectsto = 0;
             cbc1(6).onface     = 0;
           else
             cbc1(6).connectsto = ieg+1;
             cbc1(6).onface     = 5;
           end
%          Find connecting element on faces 1-4
           for k=1:4
             if ~strcmpi(cbc1(k).bc,'E  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else
               c2d_el = mesh2d.cbc(k,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + j;
               else  
                 c3d_el = j;
               end  

               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = mesh2d.cbc(k,i).onface;
             end  
           end
           cbc(:,e3d) = cbc1;
         end        % j=1:nz

       elseif mesh3d.ifzc(il)==1
%        For the layer that has been coarsened 
         for j=1:nz
           e3d=e3d+1;
           ieg = glnos(j);
%          Setting face boundary conditions
           for k=1:nfaces3
             cbc1(k).bc = mesh3d.EF{ieg}(k,:);
           end  
%          Find connecting element on face 5
           if strcmpi(cbc1(5).bc,'P  ')          % if its a periodic face
             cbc1(5).connectsto = ieg+nz-1;
             cbc1(5).onface     = 6;
           elseif ~strcmpi(cbc1(5).bc,'E  ')     % If its a boundary condition
             cbc1(5).connectsto = 0;
             cbc1(5).onface     = 0;
           else
             cbc1(5).connectsto = ieg-1;
             cbc1(5).onface     = 6;
           end  
%          Find connecting element on face 6
           if strcmpi(cbc1(6).bc,'P  ')
             cbc1(6).connectsto = ieg-nz+1;
             cbc1(6).onface     = 5;
           elseif ~strcmpi(cbc1(6).bc,'E  ')     % If its a boundary condition
             cbc1(6).connectsto = 0;
             cbc1(6).onface     = 0;
           else
             cbc1(6).connectsto = ieg+1;
             cbc1(6).onface     = 5;
           end
%          Find connecting element/onface on faces 1-3
           for k=1:3
             if ~strcmpi(cbc1(k).bc,'E  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else
               c2d_el = mesh2d.cbc(k,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + j;
               else  
                 c3d_el = j;
               end  
               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = mesh2d.cbc(k,i).onface;
             end  
           end
%          Find connecting element/onface on face 4
           k=4; 
           if ~strcmpi(cbc1(k).bc,'E  ')
             cbc1(k).connectsto = 0;
             cbc1(k).onface     = 0;
           else
             c2d_el = mesh2d.cbc(k,i).connectsto;
             if mod(j,4)==1
               c3d_el = rsum_z(c2d_el-1) + (j-1)/2 + 1;
               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = mesh2d.cbc(k,i).onface;
             elseif mod(j,4)==2
               c3d_el = rsum_z(c2d_el-1) + (j-2)/2 + 1;
               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = 6;
             elseif mod(j,4)==3
               c3d_el = rsum_z(c2d_el-1) + (j-3)/2 + 2;
               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = 5;
             else
               c3d_el = rsum_z(c2d_el-1) + (j-4)/2 + 2;
               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = mesh2d.cbc(k,i).onface;
             end   % mod(j,4)  
           end     % strcmpi(cbc1(k).bc,'E  ')
           cbc(:,e3d) = cbc1;
         end       % j=1:nz 

       elseif mesh3d.ifzc(il)==2
%        For the layer that is above the coarsened layer 
         for j=1:nz
           e3d=e3d+1;
           ieg = glnos(j);
%          Setting face boundary conditions
           for k=1:nfaces3
             cbc1(k).bc = mesh3d.EF{ieg}(k,:);
           end  
%          Find connecting element on face 1,3 and 4
           for k=[1,3,4]
             if ~strcmpi(cbc1(k).bc,'E  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else
               c2d_el = mesh2d.cbc(k,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + j;
               else  
                 c3d_el = j;
               end  
               cbc1(k).connectsto = c3d_el;
               cbc1(k).onface     = mesh2d.cbc(k,i).onface;
             end  
           end
%          Find connecting element/onface on face 2,5 and 6
           if mod(j,2)==1
%            face 2
             k=2;
             if ~strcmpi(cbc1(k).bc,'E  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else  
               c2d_el = mesh2d.cbc(k,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + (j-1)*2 + 1;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = mesh2d.cbc(k,i).onface;
               else
                 c3d_el = 0 + (j-1)*2 + 1;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = mesh2d.cbc(k,i).onface;
               end
             end  % strcmpi(cbc1(k).bc,'E  ')
%            face 5
             k=5;
             if strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = ieg+nz-1;
               cbc1(k).onface     = 6;
             elseif ~strcmpi(cbc1(k).bc,'E  ') && ~strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else  
               cbc1(k).connectsto = ieg-1;
               cbc1(k).onface     = 6;
             end
%            face 6
             k=6;
             if ~strcmpi(cbc1(k).bc,'E  ') && ~strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else  
               c2d_el = mesh2d.cbc(2,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + (j-1)*2 + 2;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = 4;
               else
                 c3d_el = 0 + (j-1)*2 + 2;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = 4;
               end
             end
           else         % mod(j,2)==0
%            face 2
             k=2;
             if ~strcmpi(cbc1(k).bc,'E  ') && ~strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else  
               c2d_el = mesh2d.cbc(k,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + j*2;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = mesh2d.cbc(k,i).onface;
               else
                 c3d_el = 0 + j*2;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = mesh2d.cbc(k,i).onface;
               end
             end  % strcmpi(cbc1(k).bc,'E  ')
%            face 5
             k=5;
             if ~strcmpi(cbc1(k).bc,'E  ') && ~strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else  
               c2d_el = mesh2d.cbc(2,i).connectsto;
               if c2d_el>1
                 c3d_el = rsum_z(c2d_el-1) + j*2 - 1;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = 4;
               else
                 c3d_el = 0 + j*2 - 1;
                 cbc1(k).connectsto = c3d_el;
                 cbc1(k).onface     = 4;
               end
             end
%            face 6
             k=6;
             if strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = ieg-nz+1;
               cbc1(k).onface     = 5;
             elseif ~strcmpi(cbc1(k).bc,'E  ') && ~strcmpi(cbc1(k).bc,'P  ')
               cbc1(k).connectsto = 0;
               cbc1(k).onface     = 0;
             else 
               cbc1(k).connectsto = ieg+1;
               cbc1(k).onface     = 5;
             end
           end   % mod(j,2)  
           cbc(:,e3d) = cbc1;
         end       % j=1:nz 

       else
         disp('Some element missing for connectivity')        
       end          % ifzc
     end            % i=lind_2d
   end              % il=1:nlayers 
      
   mesh3d.cbc = cbc;
end   % function

%----------------------------------------------------------------------



