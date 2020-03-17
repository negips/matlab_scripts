% Small test simulation.

clear
clc
close all

%addpath '../../'

SIZE           %    Polynomial order
re2            %    domain/decomposition/mapping/boundary conditions.

ifplot = 0;
gltot = 0;

display(['Time started: ' datestr(clock)])

for ii=1:nelv

     display(['Building Matrices for Element ', num2str(ii)])

     [mass nek_mass DXM1 DYM1 DXM1D DYM1D RXM1 RYM1 SXM1 SYM1 convx convy convall convxd convyd convalld gradm1xd gradm1yd intpm1d wtsvecd nek_conv lpx lpy lpall nek_lp lpbc forc x_coeff y_coeff Dx Dy w2m1 xm1 ym1 JACM1] = MESem2D4(Nx,Ny,Nxd,Nyd,El(ii).xc,El(ii).yc,ifplot);

     El(ii).mass = mass;      
     El(ii).nek_mass = nek_mass;
     El(ii).Dx = Dx;
     El(ii).Dy = Dy;
     El(ii).x_coeff = x_coeff;
     El(ii).y_coeff = y_coeff;
     El(ii).DXM1 = DXM1;
     El(ii).DYM1 = DYM1;
     El(ii).DXM1D = DXM1D;
     El(ii).DYM1D = DYM1D;
     El(ii).RXM1 = RXM1;
     El(ii).RYM1 = RYM1;
     El(ii).SXM1 = SXM1;
     El(ii).SYM1 = SYM1;
     El(ii).convx = convx;
     El(ii).convy = convy;
     El(ii).convall = convall;
     El(ii).convxd = convxd;
     El(ii).convyd = convyd;
     El(ii).convalld = convalld;
     El(ii).gradm1xd = gradm1xd;
     El(ii).gradm1yd = gradm1yd;
     El(ii).intpm1d = intpm1d;
     El(ii).wtsvecd = wtsvecd;

     El(ii).nek_conv = nek_conv;
     El(ii).lapl = lpall;
     El(ii).nek_lp = nek_lp;
     El(ii).xm1 = xm1;
     El(ii).ym1 = ym1;
     El(ii).JACM1 = JACM1;
     El(ii).w2m1 = w2m1;
%     [El(ii).lglmap gltot] = GlobalNo(ii,Nx,Ny,nelx,nely,xperiodic,yperiodic,gltot,El);

     El(ii).un = 0*xm1;
     El(ii).ulag1 = El(ii).un;
     El(ii).ulag2 = El(ii).un;
     El(ii).soln = El(ii).un;
     El(ii).b = El(ii).un;
     El(ii).bl = El(ii).un;
     El(ii).scrtch1 = El(ii).un;
     El(ii).scrtch2 = El(ii).un;

end

%    Build multiplicity
for ii = 1:nelv
     [El(ii).lglmap gltot] = GlobalNo(ii,Nx,Ny,nelx,nely,xperiodic,yperiodic,gltot,El);
end
El = vmult(El,nelv,lx1,ly1);

%    inidices/Global nos of boundary nodes
El = BdryIndex(El,nelv,lx1,ly1);

for elno=1:nelv
     for jj=0:Ny
          for ii=0:Nx
               xpt=El(elno).xm1(ii+1,jj+1);
               ypt=El(elno).ym1(ii+1,jj+1);
               El(elno).un(ii+1,jj+1) = usric(xpt,ypt);
               posx = jj*(Nx+1) + ii + 1;
               El(elno).unvec(posx,1) = El(elno).un(ii+1,jj+1);
          end
     end
end

plotgll=1;
nekchecks(El,Nx,Ny,nelx,nely,nelv,plotgll)

[bigmass bigconv bigconvd bigconvxd velvec gno nreps nn] = AssembleBig(El,Nx,Ny,nelx,nely,nelv);

clearvars mass nek_mass DXM1 DYM1 DXM1D DYM1D RXM1 RYM1 SXM1 SYM1 convx convy convall convxd convyd convalld gradm1xd intpm1d wtsvecd nek_conv lpx lpy lpall nek_lp lpbc forc x_coeff y_coeff Dx Dy w2m1 xm1 ym1 JACM1

clearvars xpt ypt ee ii jj

save(['matrices_N' num2str(Nx) '_NELV' num2str(nelv) '_CART3.mat']);
%% SOLVE

