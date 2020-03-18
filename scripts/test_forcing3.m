%    Test for creating forcing vector


clear
clc
close all

file = 'fst_spectrum.dat';

spectrum = importdata(file);
shell = spectrum.data(:,1);
kx    = spectrum.data(:,2);
ky    = spectrum.data(:,3);
kz    = spectrum.data(:,4);
u     = spectrum.data(:,5);
v     = spectrum.data(:,6);
w     = spectrum.data(:,7);
fx    = spectrum.data(:,8);
fy    = spectrum.data(:,9);
fz    = spectrum.data(:,10);

K     = spectrum.data(:,2:4);
U     = spectrum.data(:,5:7);

F     = 0*U;
ffx   = 0*u;
ffy   = 0*v;
ffz   = 0*w;

i     = 110;

nmodes=length(kx);
residuals =zeros(nmodes,1);

for i=1:nmodes

  clc
  
  vecr=U(i,:)';
  
  Re = 100000.;
  k2 = kx(i)^2 + ky(i)^2 + kz(i)^2;
  k4 = k2^2;
  
  omega = kx(i);
  
  vecr   = vecr*k2/Re;
  F(i,:) = vecr;
  ffx(i)  = vecr(1);
  ffy(i)  = vecr(2);
  ffz(i)  = vecr(3);

end 

%plot3(u,v,w,'.'); hold on
plot3(ffx,ffy,ffz,'r.'); hold on
plot3(fx,fy,fz,'k.'); hold on



