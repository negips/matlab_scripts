function [soln err iter] = me_cggo(A,x0,b,cbc,max_it,tol)

   [nx ny nz] = size(A);

   if (nz == -1)
        error(message('Single Element. Use standard GMRES./cggo'));
   end 

   [nxb nzb]  = size(b);  
   if (nzb~=nz) || (nxb~=nx)
     error(message('Dimension mismatch between A and b'));
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

   nels=nz;  

   A2 = A;
   B2 = b;
   C2 = cbc; 
   bcvec = zeros(nx,1);
   lx = x0; 
   gb = DSSUM(b,nels,1);
   Ax = zeros(nx,nels);
   relres = zeros(nels,1);   

   hdebug = figure;  

   itrs = 0;          
   while itrs<max_it
   
     itrs = itrs+1;
     for nel=1:nels
% Boundry terms.
% assuming more than 1 element.
%        bcm = C2(:,:,nel);  
%        bc = getbdry(lx,bcm,nx,nel,nels)  

        mass = A2(:,:,nel);
        b = B2(:,nel);
        b = b-bc;

        [x,flag,relres(nel)] = pcg(mass,b,tol);
        lx(:,nel) = x;
        Ax(:,nel) = mass*x;  
     end

% Do DSSM
%     lx = DSMEAN(lx,nels,1);
%     lx = DSSUM(lx,nels,1);
     [ifconv resid] = evalconv(A2,B2,lx,nx,nels);

%     lx = DSSUM(lx,nels,1);

     x=lglnodes(nx-1);
     x=x(end:-1:1);
     x = [x x+2];
     figure(hdebug);
     plot(x,lx,'-sr')
     grid on
     dbstop in me_cggo at 72 

     display(ifconv)

     if ifconv<tol;
          break;
     end

     end
     lx = DSSUM(lx,nels,1);
     soln = lx;
     err = ifconv;
     iter = itrs;

     close(hdebug)
     return
%% End of me_cggo function
%---------------------------------------------------------------------- 

function var = DSSUM(v,nels,periodic)
%% DSSUM  - as nek calls it
   
   var = v;    
   if nels==1
     if periodic
          v(1) = (v(1) + v(end));
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

function var = DSMEAN(v,nels,periodic)
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
          v(end,nel) = (v(end,nel)+v(1,nel+1))/2;

          if periodic
               v(1,1) = (v(1,1) + v(end,nels))/2;
          end
     elseif nel<nels
          v(1,nel) = v(end,nel-1);
          v(end,nel) = (v(end,nel) + v(1,nel+1))/2;

     elseif nel==nels
          v(1,nel) = v(end,nel-1);

          if periodic
               v(end,nel) = v(1,1);
          end
     end
   end         
   var = v;  
   return  
%    End of DSMEAN


function [bc] = getbdry(lx,bcm,nx,nel,nels)

     bc=zeros(nx,1);
     if nel==1
          bc(1,1) = lx(end,end);
          bc(end,1) = lx(1,nel+1);
     elseif nel==nels
          bc(1,1) = lx(end,nel-1);
          bc(end,end) = lx(1,1);
     else
          bc(1,nel) = lx(end,nel-1);
          bc(end,nel) = lx(1,nel+1);
     end
     bc = bcm*bc;
     return

function [conv r] = evalconv(A,b,x,nx,nels)

     conv =0;
     r = zeros(nx,nels);
     r1 = zeros(nx,1);
     for nel =1:nels
        Ax = A(:,:,nel)*x(:,nel); 
        if  nel==1
               r(:,1) = b(:,nel) - A(:,:,nel)*x(:,nel);
        elseif nel==nels
               r1 = b(:,nel) - A(:,:,nel)*x(:,nel);
               r(end,nel-1) = r(end,nel-1) + r1(1);
               r1(1) = 0;
               r1(end) = r1(end) + r(1,1);
               r(:,nel) = r1;
               r(1,1) = 0;
        else
               r1 = b(:,nel) - A(:,:,nel)*x(:,nel);
               r(end,nel-1) = r(end,nel-1) + r1(1);
               r1(1) = 0;
               r(:,nel) = r1;
        end  

     end

     conv=0;
     for nel=1:nels
          conv = conv + sum(r.^2);
     end
     conv = sqrt(conv);
     return

% Matrix splitting BC

