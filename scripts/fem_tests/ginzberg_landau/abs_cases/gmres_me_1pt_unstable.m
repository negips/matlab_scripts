function [soln, glerr, t_iter, flag] = gmres_me( A, x0, b, M, restrt, max_it, tol, periodic )

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


% Modified for multi dimensional load vectors
% for Multi element SEMs.     

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
   b2 = b;    
   for nel =1:nz
     mass(:,:,nel) = A2(:,:,nel);
     mass(1,:,nel) = zeros(1,nx);
     mass(1,1,nel) = 1;
     loc_x(:,nel) = x0(:,nel);
     loc_b(:,nel) = b2(:,nel); 
   end

%% Declarations
%   mass = A;  
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

%% Check initial guess   
%   gl_b = getglobal(b,nx,nz,dof);            % Global load vector
%   loc_b = b;                                % local load vectors
%   loc_x = x0;                               % local initial solution  

%   bnrm2 = norm(gl_b);
%   if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

%% Modify below
   for nel=1:nels
        b=loc_b(:,nel);     
        A=mass(:,:,nel);
        x = x0(:,nel);  
        r = b-A*x;
        lb2norm(nel) = norm(b);
       
        if lb2norm(nel)~=0
          conv(nel) =  norm(r)/lb2norm(nel);
        else
          conv(nel) = norm(r);
        end
    
        glerr(:,nel) = r;
        ifconv = zeros(1,nel);  
        if conv(nel)<=tol
          ifconv(nel) = 1;
        else
          ifconv(nel) = 0;  
        end
   end 

   if sum(ifconv)==nels
        display('Initial solution is converged');
        soln = x0;
        t_iter = 0;
        glerr = r;
        return;
   end

% if its not converged ...

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
        if (outer_itrs~=1) || (nel~=1)          % No communication in first step for first element.
          if nel==1
               b(1) = loc_x(end,end);
          else
               b(1) = loc_x(end,nel-1);
          end
        end 

%    If rhs == 0. Trivial solution for this element.
        if lb2norm(nel) == 0
          x = zeros(nx,1);
          loc_x(:,nel) = x;
          continue;
        end 

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
     loc_x(1,1) = loc_x(end,end);

%% Check global convergence including BCs
%
   ifconv = zeros(1,nels);  
   for nel=1:nels
        b=b2(:,nel);     
        A=A2(:,:,nel);
        x = loc_x(:,nel);
        r = b-A*x;
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

   lb2norm;  
   [conv tol];
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


