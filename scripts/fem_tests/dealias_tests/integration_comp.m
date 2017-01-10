% Testing 2d integration using spectral interpolation/Mass matrix
% and Full dealiasing.

addpath '../2d/'

Nx=16;
Nxd1=26; 
Nxd2=33; 

[int_fld1 spec_l1 int_fld2 spec_l2] = paul_int(Nx,Nxd1,Nxd2);

dealias_fld = dealias_int(Nx,Nxd1,Nxd2); 

[x wx p] = lglnodes(Nx);
x    =    x(end:-1:1);
wx   =    wx(end:-1:1);

figure
plot(x,int_fld1);
hold on
plot(x,dealias_fld, 'k')

figure
plot(x,int_fld1./dealias_fld,'k')

