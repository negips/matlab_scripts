function [soln err iter] = cggo_me1(A,x0,b,cbc,max_it,tol)

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
   lx = x0;
   lx2 = 0*x0;   
   gb = DSSUM(b,nels,1);
   Ax = zeros(nx,1);

   bcm1 = zeros(nx,nx,nels);
   bcm2 = zeros(nx,nx,nels);
   bcvec1 = zeros(nx,1);
   bcvec2 = zeros(nx,1);     
   lvec = zeros(nx,1);  

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

   hdebug = figure;  

   itrs = 0;          
   while itrs<max_it
   
     itrs = itrs+1;
     for nel=1:nels

% Add force terms due to split matrix
% assuming more than 1 element.
%        bcm = C2(:,:,nel);  
%        bc = getbdry(lx,bcm,nx,nel,nels)  

        [bcvec1 bcvec2] = getforce(bcm1,bcm2,nel,nels,lx);
%        lvec = getloadvec(B2,nel,nels);   

        mass = A2(:,:,nel);
        b = gb(:,nel);
%        b = b-bc;
        b = b - bcvec1 - bcvec2;  

        [x,flag,relres] = pcg(mass,b,1e-10);
        lx2(:,nel) = x;
     end

%    Check convergence
     lx = lx2;
     for nel=1:nels
          Ax = A2(:,:,nel)*lx(:,nel);
          r = B2(:,nel) - Ax;
          rr = norm(r)/norm(B2(:,nel));
          if rr<tol
               ifconv(nel) = 1;
          else
               ifconv(nel) = 0;
          end
     end

%     [ifconv resid] = evalconv(A2,B2,lx,nx,nels);

     x=lglnodes(nx-1);
     x=x(end:-1:1);
     x = [x x+2];
     figure(hdebug);
     plot(x,lx,'-sr')
     hold on
     grid on
%     dbstop in cggo_me at 94 

     display(ifconv)

     if sum(ifconv)==nels;
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
               r2 = -A(:,:,nel)*x(:,nel);
               r1(1) = r2(1);
               r(end,nel-1) = r(end,nel-1) + r1(1);
               r1(1) = 0;
               r1(end) = r1(end) + r(1,1);
               r(:,nel) = r1;
               r(1,1) = 0;
        else
               r1 = b(:,nel) - A(:,:,nel)*x(:,nel);
               r2 = -A(:,:,nel)*x(:,nel);
               r1(1) = r2(1);
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

function lv = getloadvec(b,nel,nels)

     [nx nz] = size(b);     
     lv = zeros(nx,1);
     if nel==1
          lv(1) = b(end,end);
          lv(end) = b(1,2);
     elseif nel==nels
          lv(1) = b(end,nel-1);
          lv(end) = b(1,nel-1);
     else
          lv(1) = b(end,nel-1);
          lv(end) = b(1,nel+1);
     end
     return








