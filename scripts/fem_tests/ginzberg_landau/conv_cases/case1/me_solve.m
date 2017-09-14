function [soln, glerr, t_iter, flag] = me_solve( A,x0,b, max_it, tol)

%  -- Iterative template routine --

%  multi dimensional load vectors/Mass matricies.
%  for Multi element SEMs.     

   [nx ny nz] = size(A);

   if (nz == -1)
        error(message('Single Element. Use standard GMRES.'));
   end 

   [nxb nzb]  = size(b);  
   if (nzb~=nz) || (nxb~=nx)
     error(message('Dimension mismatch between A and b'));
   end

%    Check initial guess dimensions
   if ~isempty(x0)
        [nxx nzx]  = size(x0);  
        if (nzx~=nz) || (nxx~=nx)
          error(message('Dimension mismatch between A and x0'));
        end
   else
        x0 = zeros(nx,nz);
   end    
%    Check zero load vector.

   isloadzero = max(max(abs(b)));
   if isloadzero==0
     display('Trivial solutions for {b} = 0');
     soln = b;
     t_iter=0;
     glerr = 0;
     flag = 1;
     return
   end  

%% Initialization
% Since the first point acts as a boundry...
   A2 = A;
   B2 = b;
   mass = zeros(nx,ny,nz);
   lx = x0;
   lb = B2;        

%% Declarations
   mass = A2;  
   nx = nx;
   ny = ny;    
   nels = nz;
   dof = (nx-1)*nz+1;

   t_iter = 0;
   flag = 0;
   max_outer_itrs = max_it;  

   conv = zeros(1,nels);
   soln = zeros(nx,nels);    
   bcvec = zeros(nx,1);
   bcm1=zeros(1,nx);
   bcm2=zeros(1,nx);    
  
%% Check initial guess   

% if its not converged ...

   for nel=1:nels
     if nel==1
          bcm1(1,:,nel) = A2(end,:,nels);
          bcm2(end,:,nel) = A2(1,:,2);
     elseif nel==nels
          bcm1(1,:,nel) = A2(end,:,nel-1);
          bcm2(end,:,nel) = A2(1,:,1);
     else
          bcm1(1,:,nel) = A2(end,:,nel-1);
          bcm2(end,:,nel) = A2(1,:,nel+1);
     end
   end   


%% SOLVE
   outer_itrs = 0;
   tol = 1e-9;  
  
   while (outer_itrs<max_outer_itrs)

        outer_itrs = outer_itrs+1;  
     
        for nel=1:nels

             b=lb(:,nel);
             A=mass(:,:,nel);  

%         Communication with neighboring elements
             [bcvec1 bcvec2] = getforce(bcm1,bcm2,nel,nels,lx);
             b(1) = b(1) - bcvec1;
             b(end) = b(end)-bcvec2;         

             [x,cflag,relres] = pcg(A,b,tol);

             lx(:,nel) = x;   

        end           % nel=1:nels

        [glx lx]=gllosolution(lx); 

%        plot(glx)  
%dbstop in me_solve at 112

   end           % while ~convergence
   soln = lx;
   glerr = 0;
   t_iter = outer_itrs;
   flag = 0; 
   return              
% END of me_solve 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix splitting as forcing 
function [bcvec1 bcvec2] = getforce(bcm1,bcm2,nel,nels,lx)

     if nel==1
          bcvec1 = bcm1(:,:,nel)*lx(:,end);
          bcvec2 = bcm2(:,:,nel)*lx(:,2);
     elseif nel==nels
          bcvec1 = bcm1(:,:,nel)*lx(:,nel-1);
          bcvec2 = bcm2(:,:,nel)*lx(:,1);
     else
          bcvec1 = bcm1(:,:,nel)*lx(:,nel-1);
          bcvec2 = bcm2(:,:,nel)*lx(:,nel+1);
     end
     return

function [gl_soln lo_soln] = gllosolution(lx)

     [nx nz] = size(lx);
      dof = nz*(nx-1)+1;
      xg=zeros(dof,1);
      xl=0*lx;

      for nel=1:nz
         if nel==1
               xg(1:nx)=lx(:,nel);
               count=nx;
               xl=lx(:,nel);
          else
               xg(count)=(xg(count)+lx(1,nel))/2;
               xg(count+1:count+nx-1)=lx(2:end,nel);
               count=count+nx-1;
               xl(:,nel) = lx(:,nel);
               xl(1,nel) = (xl(1,nel)+lx(end,nel-1))/2;
               xl(end,nel-1)=xl(1,nel);
          end
      end
      xg(1)=(xg(1)+xg(end))/2;
      xg(end) = xg(1);
      gl_soln=xg;

      xl(end,nz)=(xl(end,nz)+xl(1,1))/2;
      xl(1,1)=xl(end,nz);
      lo_soln=xl;

      return





