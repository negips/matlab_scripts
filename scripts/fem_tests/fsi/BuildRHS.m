function [rhs El]= BuildRHS(El,Nx,Ny,gno,nn,nelv)

     lbasis = (Nx+1)*(Ny+1);
     convecting_fld = El(1).unvec*0 + 1;
     rhs = zeros(nn,1);

     for elno=1:nelv

          c1 = El(elno).intpm1d*convecting_fld;                                 % Convecting field
          d1 = (El(elno).gradm1xd +0*El(elno).gradm1yd)*El(elno).unvec;           % derivate
          rhsvec = zeros(lbasis,1);
    
          for ii=0:lbasis-1
               test_fn = zeros(lbasis,1);
               test_fn(ii+1) = 1;

               t1 = El(elno).intpm1d*test_fn;                              % test function
               mul = El(elno).wtsvecd.*c1.*d1.*t1;                         % weight*Product 
               summ = sum(mul);

               posx = mod(ii,Nx+1)+1;
               posy = floor(ii/(Nx+1))+1;
               El(elno).rhs(posx,posy) = summ;

               gl_no = El(elno).lglmap(posx,posy);
               ind = find(gno==gl_no);
               rhs(ind) = rhs(ind) + summ; 
          end

     end

%     dbstop in BuildRHS at 33
     return

