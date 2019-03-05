function [soln, glerr, iter, flag] = gmres_me( A, x0, b, M, restrt, max_it, tol, periodic )

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

%    Check initial gues dimensions
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
     iter=0;
     glerr = zeros(nx,nz);
     flag = 1;
     return
   end  



%%  Initialization
   mass = A;  
   nels = nz;
   dof = (nx-1)*nz+1;

   iter = 0;
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

%% Check initial guess   
%   gl_b = getglobal(b,nx,nz,dof);            % Global load vector
   loc_b = b;                                % local load vectors
   loc_x = x0;                               % local initial solution  

%   bnrm2 = norm(gl_b);
%   if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
   Ax = zeros(nx,nels);  
   l_r = zeros(nx,nels);  

   for nel=1:nels
        b=loc_b(:,nel);     
        A=mass(:,:,nel);
        x = loc_x(:,nel);  
        r = b-A*x;
        glerr(:,nel) = r;           
%        l_r(:,nel) = b-A*loc_x(:,nel);
        lb2norm(nel) = norm(b);
        conv(nel) =  norm(r)/lb2norm(nel);
        if conv(nel)<=tol
          conv(nel) = 1;
        end     
   end 
   
   if sum(conv)==nels
        display('Initial solution is converged');
        soln = x0;
        iter = 0;
        glerr = r;
        return;
   end
       
%   gl_r = getglobal(l_r,nx,nz,dof);          % global residual
%   glerr = norm( gl_r )/bnrm2;

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
        if lb2norm(nel) == 0
          x = b;
          continue;
        end 

        A=mass(:,:,nel);
        x = loc_x(:,nel);  
        r = b-A*x;
%        r = M \ ( b-A*x );
           
        lerr = norm(r)/lb2norm(nel);
        if ( lerr <= tol )
          conv(nel) = 1; 
          continue;                          % Converged. Go to next element.
        end

        for iter = 1:max_it                  % begin outer iterations for element.
 
          r = b-A*x;
          V(:,1) = r / norm( r );            % normalized residual
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
                inner_conv_n = i; 

                if ( err <= tol )                                % update approximation
                      y = H(1:i,1:i) \ s(1:i);                % and exit inner iterations
                      x = x + V(:,1:i)*y;
                      loc_x(:,nel) = x;
                      conv(nel) = 1;      
                   break;
                end

           end           % i=1:m        % inner iterations
          
           if (err<=tol)
               break;
           else
                y = H(1:m,1:m) \ s(1:m);
                x = x + V(:,1:m)*y;                    % update approximation for restart
                loc_x(:,nel) = x;
           end 

        end                                       % iter = 1:max_it for each element.

%        loc_x = EleConvOp(loc_x,nel,nels,periodic);              % for now no passing between solves. 

   end           % nel=1:nels

   loc_x = DSSUM(loc_x,nels,periodic);

   for nel=1:nels
        b=loc_b(:,nel);     
        A=mass(:,:,nel);
        x = loc_x(:,nel);  
        r = b-A*x;         
        conv(nel) =  norm(r)/lb2norm(nel);
        glerr(:,nel) = r;
        if conv(nel)<=tol
          conv(nel) = 1;
        end     
   end 

   soln = loc_x;
   iter = outer_itrs;
   flag=0;   
   
   if sum(conv)==nels
        display('Solution is converged');
        flag=1;   
        break;
   end

   end           % while ~convergence    
% END of gmres.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function var = DSSUM(v,nels,periodic)
%% DSSUM  - as nek calls it
   
   var = v;    
   if nels==1
     if periodic
          v(1) = (v(1) + v(end))/2;
          v(end) = v(1);
          var = v;
          return;
     else
          return;
     end 
   end   

   for nel = 1:nels
     if nel==1
          v(end,nel) = v(end,nel)+v(1,nel+1);

          if periodic
               v(1,1) = v(1,1) + v(end,nels);
          end
     elseif nel<nels
          v(1,nel) = v(end,nel-1);
          v(end,nel) = v(end,nel) + v(1,nel+1);

     elseif nel==nels
          v(1,nel) = v(end,nel-1);

          if periodic
               v(end,nel) = v(1,1);
          end
     end
   end         
   var = v;  
   return  
%    End of DSSUM


function var = EleConvOp(lvar,nel,nels,periodic)
%%   Communication between elements
     var = lvar;
     if nels==1
          if periodic
               lvar(1) = (lvar(1)+lvar(end))/2;
               lvar(end) = lvar(1);
               var = lvar;
               return;
          end
     end

     if nel~=nels
          lvar(1,nel+1) = lvar(1,nel+1) + lvar(end,nel);
          var = lvar;
     else
          if periodic
               lvar(1,1) = (lvar(1,1)+lvar(end,nels))/2;
          end
          var = lvar;
     end
               
               


function glvar = getglobal(lvar,nx,nz,dof)
     
   glvar = zeros(dof,1);           % global variable 
   for i=1:(nz-1)
        glvar((i-1)*(nx-1)+1:(i-1)*(nx-1)+nx) = glvar((i-1)*(nx-1)+1:(i-1)*(nx-1)+nx)+lvar(:,i);
   end
   glvar((nz-1)*(nx-1)+1:(nz-1)*(nx-1)+nx)= glvar((nz-1)*(nx-1)+1:(nz-1)*(nx-1)+nx) + lvar(:,nz);


