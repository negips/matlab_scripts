function [soln, glerr, t_iter, flag] = gmres_me( A, x0, b, bc, M, restrt, max_it, tol, periodic )

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag] = gmres( A, x, b, M, restrt, max_it, tol )
%
% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   A        REAL nonsymmetric positive definite matrix
%         x0       REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         restrt   INTEGER number of iterations between restarts
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  soln     REAL solution vector
%         err      REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it


% Modified for multi dimensional load vectors/Mass matricies.
% for Multi element SEMs.     

% Using Boundary terms for the convection equation.

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
     glerr = zeros(nx,nz);
     flag = 1;
     return
   end  

%% Initialization
% Since the first point acts as a boundry...
   A2 = A;
   B2 = b;
   CBC2 = bc;
   mass = zeros(nx,ny,nz);
   loc_bc = zeros(nx,ny,nz);
   loc_x = x0;
   loc_b = B2;        

%% Declarations
   mass = A2;  
   nx = nx;
   ny = ny;    
   nels = nz;
   dof = (nx-1)*nz+1;

   t_iter = 0;
   flag = 0;
   max_outer_itrs = 20;  

   m = restrt;  
   V = zeros(nx,m+1);
   H = zeros(m+1,m);
   cs = zeros(m,1);
   sn = zeros(m,1);
   e1 = zeros(m,1);
   e1(1) = 1;
   s = zeros(m,1);
   conv = zeros(1,nels);
   soln = zeros(nx,nels);    
   bcvec = zeros(nx,1);
  
%% Check initial guess   
%   gl_b = getglobal(b,nx,nz,dof);            % Global load vector
%   loc_b = b;                                % local load vectors
%   loc_x = x0;                               % local initial solution  



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
  
   while (outer_itrs<max_outer_itrs)
   outer_itrs = outer_itrs+1;  
     
   for nel=1:nels

        V = zeros(nx,m+1);
        H = zeros(m+1,m);
        cs = zeros(m,1);
        sn = zeros(m,1);

        b=loc_b(:,nel);

%    Build boundary conditions
       [bcvec1 bcvec2] = getforce(bcm1,bcm2,nel,nels,lx)

        b = b - bcvec1 - bcvec2;                            % Now this has information of boundaries

        A=mass(:,:,nel);
        x = loc_x(:,nel);  
        r = b-A*x;
%        r = M \ ( b-A*x );

%    If initial solution is good enough.           
        lerr = norm(r)/lb2norm(nel);
        if ( lerr <= tol )
          conv(nel) = norm(r); 
          loc_x(:,nel) = x;
          continue;                          % Converged. Go to next element.
        end

        for iter = 1:max_it                  % begin outer iterations for element.
 
          r = b-A*x;
          V(:,1) = r / norm( r );            % normalized krylov vector?
          s = norm( r )*e1;

     %    Inner iterations 
          for i = 1:m                            % construct 'm' orthonormal vector
          %	 w = M \ (A*V(:,i));                    % basis using Gram-Schmidt
                w = A*V(:,i);                    % next Krylov space vector (non-orthogonal)

                for k = 1:i                            % Orthonormalization procedure
                     H(k,i)= w'*V(:,k);          % Hessenberg matrix
                     w = w - H(k,i)*V(:,k);      % residual after projection to kth vector space
                end
                H(i+1,i) = norm( w );
                V(:,i+1) = w / H(i+1,i);         % normalized orthogonal vector (q)

                for k = 1:i-1                               % apply Givens rotation
                      temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
                      H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
                      H(k,i)   = temp;
                end

                [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) );  % form i-th rotation matrix
                temp   = cs(i)*s(i);                                 % approximate residual norm
                s(i+1) = -sn(i)*s(i);
                s(i)   = temp;
                H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
                H(i+1,i) = 0.0;

                err  = abs(s(i+1)) / lb2norm(nel);

                if ( err <= tol )                                % update approximation
                      y = H(1:i,1:i) \ s(1:i);                % and exit inner iterations
                      x = x + V(:,1:i)*y;
                      break;
                end

           end           % i=1:m        % inner iterations
          
           if (err<=tol)
               break;
           end

           y = H(1:m,1:m) \ s(1:m);
           x = x + V(:,1:m)*y;                    % update approximation for restart
           
           r = b-A*x;
           s(i+1) = norm(r);
           err = s(i+1)/lb2norm(nel);
           if err<=tol
               break;
           end 
          
        end                                       % iter = 1:max_it for each element.

        loc_x(:,nel) = x;

   end           % nel=1:nels

%    Final communication between first and last element.
%     loc_x(1,1) = loc_x(end,end);

%% Check global convergence including BCs
%
   ifconv = zeros(1,nels); 
   for nel=1:nels
        b=b2(:,nel);
        A=A2(:,:,nel);
        x = loc_x(:,nel);
        bc = loc_bc(:,:,nel)*x;  
        r = b-bc-A*x;
        glerr(:,nel) = r;
        soln(:,nel) = x;    
      
        if lb2norm(nel)~=0            
          conv(nel) =  norm(r)/lb2norm(nel);
        else
          conv(nel) = norm(r);
        end
          
        if conv(nel)<=tol
          ifconv(nel) = 1;
        else
          ifconv(nel) = 0;  
        end     
   end 

   t_iter=outer_itrs;
   flag=0; 
  
   if sum(ifconv)==nels
        display('Solution is converged');
        flag=1;   
        break;
   end

   end           % while ~convergence
% END of gmres.m
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
