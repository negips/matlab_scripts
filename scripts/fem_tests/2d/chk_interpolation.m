% Verifying nek's interpolation

clear
clc

%% Over-integration matrices using higher order GL points (Gauss-Legendre)
%-------------------------------------------------- 
Nx=2;
Ny=Nx;
Nxd=4;
Nyd=Nxd;

if Nx>0
     [x wx p]= lglnodes(Nx);
     x =x(end:-1:1);
     wx =wx(end:-1:1);
else
     x=0;
     wx=1;
end

if Ny>0
     [y wy p]= lglnodes(Ny);
     y =y(end:-1:1);
     wy = wy(end:-1:1);
else
     y=0;
     wy=1;
end


if Nxd>1
  [xGL wxGL]= lgwt(Nxd+1,-1,1);
  xGL = xGL(end:-1:1);
  wxGL= wxGL(end:-1:1);
else
  xGL=[0];
  wxGL=[1];
end

if Nyd>1
  [yGL wyGL]= lgwt(Nyd+1,-1,1);
  yGL = yGL(end:-1:1);
  wyGL= wyGL(end:-1:1);
else
  yGL=[0];
  wyGL=[1];
end

%% Interpolate to Gauss Legendre points.
%% Work in progress

%% GLL spectral to Gauss Legendre nodal
%% In the x-direction
pht = zeros(Nxd+1,Nx+1);
for j = 1:Nxd+1
  z=xGL(j);
  Lj = legendrePoly(Nx,z);
  pht(j,:) = transpose(Lj);
end
Nx_spec2Gauss = pht;

%% In the y-direction
pht = zeros(Nyd+1,Ny+1);
for j = 1:Nyd+1
  z=yGL(j);
  Lj = legendrePoly(Ny,z);
  pht(j,:) = transpose(Lj);
end
Ny_spec2Gauss = pht;
%-------------------------------------------------- 
pht = zeros(Nxd+1,Nxd+1);
for j = 1:Nxd+1
  z=xGL(j);
  Lj = legendrePoly(Nxd,z);
  pht(j,:) = transpose(Lj);
end
Nxd_spec2Gauss = pht;
Nxd_Gauss2spec = inv(pht);

%% In the y-direction
pht = zeros(Nyd+1,Nyd+1);
for j = 1:Nyd+1
  z=yGL(j);
  Lj = legendrePoly(Nyd,z);
  pht(j,:) = transpose(Lj);
end
Nyd_spec2Gauss = pht;
Nyd_Gauss2spec = inv(pht);
%-------------------------------------------------- 
% filter
filterx=eye(Nxd+1);
if (Nxd>Nx)
  for ii=Nx+2:Nxd+1
    filterx(ii,ii)=0;
  end
end  
filtery=eye(Nyd+1);
if (Nyd>Ny)
  for ii=Ny+2:Nyd+1
    filtery(ii,ii)=0;
  end
end 
%--------------------------------------------------
% Gauss spectral to GLL
pht = zeros(Nx+1,Nxd+1);
for j = 1:Nx+1
  z=x(j);
  Lj = legendrePoly(Nxd,z);
  pht(j,:) = transpose(Lj);
end
Nxd_spec2GLL = pht;

%% In the y-direction
pht = zeros(Ny+1,Nyd+1);
for j = 1:Ny+1
  z=y(j);
  Lj = legendrePoly(Nyd,z);
  pht(j,:) = transpose(Lj);
end
Nyd_spec2GLL = pht;

%--------------------------------------------------  
truncated_interpolationx = Nxd_spec2GLL*filterx*Nxd_Gauss2spec;
truncated_interpolationy = Nyd_spec2GLL*filtery*Nyd_Gauss2spec;

filtered_Gauss2GLL = kron(truncated_interpolationy,truncated_interpolationx);
%% Over-integration matrices using higher order GL points (Gauss-Legendre)
%-------------------------------------------------- 
%% Interpolate to Gauss Legendre points.
%% Work in progress

%% GLL spectral to Gauss Legendre nodal
%% In the x-direction
pht = zeros(Nxd+1,Nx+1);
for j = 1:Nxd+1
  z=xGL(j);
  Lj = legendrePoly(Nx,z);
  pht(j,:) = transpose(Lj);
end
Nx_spec2Gauss = pht;

%% In the y-direction
pht = zeros(Nyd+1,Ny+1);
for j = 1:Nyd+1
  z=yGL(j);
  Lj = legendrePoly(Ny,z);
  pht(j,:) = transpose(Lj);
end
Ny_spec2Gauss = pht;
%-------------------------------------------------- 
pht = zeros(Nxd+1,Nxd+1);
for j = 1:Nxd+1
  z=xGL(j);
  Lj = legendrePoly(Nxd,z);
  pht(j,:) = transpose(Lj);
end
Nxd_spec2Gauss = pht;
Nxd_Gauss2spec = inv(pht);

%% In the y-direction
pht = zeros(Nyd+1,Nyd+1);
for j = 1:Nyd+1
  z=yGL(j);
  Lj = legendrePoly(Nyd,z);
  pht(j,:) = transpose(Lj);
end
Nyd_spec2Gauss = pht;
Nyd_Gauss2spec = inv(pht);
%-------------------------------------------------- 
% filter
filterx=eye(Nxd+1);
if (Nxd>Nx)
  for ii=Nx+2:Nxd+1
    filterx(ii,ii)=0;
  end
end  
filtery=eye(Nyd+1);
if (Nyd>Ny)
  for ii=Ny+2:Nyd+1
    filtery(ii,ii)=0;
  end
end 
%--------------------------------------------------
% Gauss spectral to GLL
pht = zeros(Nx+1,Nxd+1);
for j = 1:Nx+1
  z=x(j);
  Lj = legendrePoly(Nxd,z);
  pht(j,:) = transpose(Lj);
end
Nxd_spec2GLL = pht;

%% In the y-direction
pht = zeros(Ny+1,Nyd+1);
for j = 1:Ny+1
  z=y(j);
  Lj = legendrePoly(Nyd,z);
  pht(j,:) = transpose(Lj);
end
Nyd_spec2GLL = pht;

%--------------------------------------------------  
truncated_interpolationx = Nxd_spec2GLL*filterx*Nxd_Gauss2spec;
truncated_interpolationy = Nyd_spec2GLL*filtery*Nyd_Gauss2spec;

filtered_Gauss2GLL = kron(truncated_interpolationy,truncated_interpolationx);
%--------------------------------------------------







