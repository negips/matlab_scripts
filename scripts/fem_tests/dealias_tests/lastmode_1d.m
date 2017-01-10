%    Testing integration with mass matricies

clear
clc
close all

addpath '../2d/'

Nx_vec = 4:15;

legs= {'N', 'N+1', 'N+2', 'N-1' 'N-2'};
for jj=1:length(Nx_vec)

     Nx = Nx_vec(jj);
     Nxd1= Nx+1;       % dealiasing points.
     Nxd2= Nx+2;
     Nxd3= Nx-1;
     Nxd4= Nx-2;
     ifplot=1;

     leg_norm = 2/(2*Nx+1);

     if Nx>0
          [x wx p]= lglnodes(Nx);
          x =x(end:-1:1);
          wx =wx(end:-1:1);
     end

     if Nxd1>0
          [xd wxd p]= lglnodes(Nxd1);
          xd1 =xd(end:-1:1);
          wxd1 =wxd(end:-1:1);
     end

     if Nxd2>0
          [xd wxd p]= lglnodes(Nxd2);
          xd2 =xd(end:-1:1);
          wxd2 =wxd(end:-1:1);
     end

     if Nxd3>0
          [xd wxd p]= lglnodes(Nxd3);
          xd3 =xd(end:-1:1);
          wxd3 =wxd(end:-1:1);
     end

     if Nxd4>0
          [xd wxd p]= lglnodes(Nxd4);
          xd4 =xd(end:-1:1);
          wxd4 =wxd(end:-1:1);
     end

     % Just taking legendre mode of order Nx
     integral =0;
     for i=1:Nx+1
          xpt=x(i);
          wts=wx(i);

          Pn = legendrePoly(Nx,xpt);
          integral = integral + wts*(Pn(end)^2);
     end
     mass(jj) = integral/leg_norm;

     % dealiased integral 1
     integrald=0;
     for i=1:Nxd1+1
          xpt=xd1(i);
          wts=wxd1(i);

          Pn = legendrePoly(Nx,xpt);
          integrald = integrald + wts*(Pn(end)^2);

     end
     deal1(jj) = integrald/leg_norm;

     % dealiased integral 2 
     integrald=0;
     for i=1:Nxd2+1
          xpt=xd2(i);
          wts=wxd2(i);

          Pn = legendrePoly(Nx,xpt);
          integrald = integrald + wts*(Pn(end)^2);
     end
     deal2(jj) = integrald/leg_norm;

     % dealiased integral 3 
     integrald=0;
     for i=1:Nxd3+1
          xpt=xd3(i);
          wts=wxd3(i);

          Pn = legendrePoly(Nx,xpt);
          integrald = integrald + wts*(Pn(end)^2);
     end
     deal3(jj) = integrald/leg_norm;

     % dealiased integral 4 
     integrald=0;
     for i=1:Nxd4+1
          xpt=xd4(i);
          wts=wxd4(i);

          Pn = legendrePoly(Nx,xpt);
          integrald = integrald + wts*(Pn(end)^2);
     end
     deal4(jj) = integrald/leg_norm;

end

plot(Nx_vec,mass,'b');
hold on
plot(Nx_vec,deal1,'-r');
plot(Nx_vec,deal2,'-dk');
plot(Nx_vec,deal3,'-m');
plot(Nx_vec,deal4,'-.g');

legend(legs);

figure
plot(Nx_vec, mass./deal1, 'o-b')





