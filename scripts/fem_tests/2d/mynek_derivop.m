% Small test simulation.

clear
clc
close all
%
%addpath '../../'

nekinit4
%load 'matrices_N11_def'
%nbasis = length(velvec);
nbasis = nn;

un    = velvec;
ulag1 = zeros(size(velvec));
ulag2 = zeros(size(velvec));
b = zeros(size(velvec));
blag1 = zeros(size(velvec));
blag2 = zeros(size(velvec));
soln = zeros(size(velvec));
minv = zeros(size(velvec));

hini=figure;
ifsolnplot=1;
Elout = UpdSoln(El,velvec,gno,nelv,lx1,ly1,hini,ifsolnplot);


%% Time integration coefficients
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];
%%

%% Rea file parameters.
deltat = 0.01;
istep = 0;
nsteps = 8000;
iostep = 10;

chi = -0.0;
re=1e5;
nu=1/re;
nu = 0;   % for now;

%% SOLVE
hsoln=figure;
time = 0;
iocount=0;
tic
for i = 0:nsteps

     if istep>0
          disp(['Step, Time, Relative residual,iter: ', num2str(istep), ', ' num2str(time), ', ' ...
          num2str(relres) ', ' num2str(pcgiter)]);     
     end

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
     [rhs El] = BuildRHS(El,Nx,Ny,gno,nbasis,nelv); 

     A = bdfk(1)*bigmass;

     b11 = bdfk(2)*un; 
     b12 = bdfk(3)*ulag1; 
     b13 = bdfk(4)*ulag2;
     b1 = bigmass*(b11+b12+b13);
%%     b1 = b11+b12+b13;
%
%     cop = bigconvxd*un;
     cop = rhs;
     b21 = extk(2)*cop;
     b22 = extk(3)*blag1;
     b23 = extk(4)*blag2;
     b2 = (b21 + b22 + b23);
%%     b2 = b21+b22+b23;

     b = -b1 -b2*deltat;

%     [bigconvxd*un cop]

%     [El ifconv] = SolveSem2D(El,lx1,ly1,nelv,max_it,tol);
%     [soln cflag iter relres] = gmres(A,b,restrt,tol,max_it);
%     cflag

     tol=1e-8;
     max_it=100;
     [soln cflag relres pcgiter] = pcg(A,b,tol,max_it);

%    Update and plot solution

     if mod(istep,iostep)==0
          ifsolnplot=1;
     else
          ifsolnplot=0;
     end
     
     El = UpdSoln(El,soln,gno,nelv,lx1,ly1,hsoln,ifsolnplot);

     if mod(istep,iostep)==0
          zlim([-1.2 1.2]);
          azim=-20;
          elev=60;
          view([azim elev]);
          if istep==0
               pause(2)
          else
               pause(0.001)
          end
          iocount=iocount+1;
          hmov(iocount) = getframe(hsoln);
     end

%    Update lag veclocities
     ulag2 = ulag1;
     ulag1 = un;
     un = soln;

%    Update lag terms for forcing/convection
     blag2 = blag1;
     blag1 = cop;

end
tt = toc;
disp(['Total Solve Time: ' num2str(tt)]);
%movie2avi(hmov,'WaveMovie_def3.avi');

%-------------------------------------------------- 
%function SetVel(El,soln,gno)
%
%     nelv = size(El);
%     [lx1 ly1] = size(El(1).xm1);
%     [~ nn] = size(soln);
%     for elno=1:nelv
%          for jj=1:ly1
%          for ii=1:lx1
%               pos = find(gno == El(elno).lglmap(ii,jj));
%               El(elno).soln(ii,jj) = soln(pos);
%          end
%          end
%     end
%     return

