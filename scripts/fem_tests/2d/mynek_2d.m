% Small test simulation.

%clear
%clc
%close all
%
%addpath '../../'

nekinit

%% Solve
%% Time integration coefficients
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];
%%

deltat = 0.01;
istep = 0;
nsteps = 5;

chi = -0.0;
re=1e5;
nu=1/re;
nu = 0;   % for now;

%% Some more matrix building
for elno=1:nelv
     for jj=1:ly1
          El(elno).dymass(:,:,jj) = diag(El(elno).nek_mass(:,jj)) - nu*El(elno).nek_lp(:,:,jj);
     end
     El(elno).soln = El(elno).un;       % Just as initialization at t=0
end

time = 0;

%% SOLVE
hsoln=figure;

for i = 1:nsteps

     istep = istep+1;
     time = time+deltat; 

     % set time integration operators
     if istep == 1
          bdfk = bdf1;
          extk = ex0;
     elseif istep == 2
          bdfk = bdf2;
          extk = ex1;
     else
          bdfk = bdf3;
          extk = ex2;
     end


%         Build right hand side
     for elno=1:nelv
          for jj=1:ly1          
               El(elno).A(:,:,jj) = bdfk(1)*El(elno).dymass(:,:,jj);

               b11 = bdfk(2)*El(elno).un(:,jj); 
               b12 = bdfk(3)*El(elno).ulag1(:,jj); 
               b13 = bdfk(4)*El(elno).ulag2(:,jj);
               b1 = El(elno).dymass(:,:,jj)*(b11+b12+b13);

               b21 = extk(2)*El(elno).un(:,jj);
               b22 = extk(3)*El(elno).ulag1(:,jj);
               b23 = extk(4)*El(elno).ulag2(:,jj);
               b2 = deltat*El(elno).nek_conv(:,:,jj)*(b21 + b22 + b23);

               El(elno).bl(:,jj) = -b1 - b2;
               El(elno).b(:,jj) = -b1 -b2;
          end
     end

     El = DSSM2D_b(El,nelv);

     max_it = 50;
     tol = 1e-9;

%     if istep==2
%          break
%     end

     [El ifconv] = SolveSem2D(El,lx1,ly1,nelv,max_it,tol);

     for elno=1:nelv
          figure(hsoln)
          surf(El(elno).xm1,El(elno).ym1,El(elno).soln);
          colorbar;
          hold on
     end
     hold off

%    Update lag veclocities
     for elno=1:nelv
          El(elno).ulag2 = El(elno).ulag1;
          El(elno).ulag1 = El(elno).un;
          El(elno).un = El(elno).soln;
     end

     pause(2)
%    Plot the solution

end

%-------------------------------------------------- 
%function SolnDisplay(El,nelv,h)
%
%     for elno=1:nelv
%          figure(h)
%          hold on
%          surf(El(elno).xm1,El(elno).ym1,El(elno).soln);
%          colorbar;
%     end
%     return

