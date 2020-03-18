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

K     = spectrum.data(:,2:4);
U     = spectrum.data(:,5:7);

F     = 0*U;

i     = 110;

nmodes=length(kx);
residuals =zeros(nmodes,1);

for i=1:1 %nmodes

  clc

  M = zeros(3,3);
  
  M(1,1) = -(ky(i)^2 + kz(i)^2);
  M(2,2) = -(kx(i)^2 + kz(i)^2);
  M(3,3) = -(kx(i)^2 + ky(i)^2);
  
  for j=1:3
    for k=1:3
      if (j ~= k)
        M(j,k) = K(i,j)*K(i,k);
      end  
    end
  end
  
  vecr=U(i,:)';
  veci=U(i,:)';
  
  Re = 1.;
  k2 = kx(i)^2 + ky(i)^2 + kz(i)^2;
  k4 = k2^2;
  
  omega = kx(i);
  
  vecr = vecr*k4/Re;

% Do QR factorization of M

  c = M(2,1);
  s = M(3,1);
  r = sqrt(c^2 + s^2);

  GV = eye(3,3);
  GV(2,2) = c/r;
  GV(2,3) = s/r;
  GV(3,2) = -s/r;
  GV(3,3) = c/r;

  Usvd = GV';

  M2 = GV*M;

  c = M2(1,1);
  s = M2(2,1);
  r = sqrt(c^2 + s^2);

  GV = eye(3,3);
  GV(1,1) = c/r;
  GV(1,2) = s/r;
  GV(2,1) = -s/r;
  GV(2,2) = c/r;

  M3 = GV*M2;

  Usvd  = Usvd*GV';
%
  c = M3(1,2);
  s = M3(1,3);
  r = sqrt(c^2 + s^2);

  GV = eye(3,3);
  GV(2,3) = -s/r;
  GV(2,2) = c/r;
  GV(3,2) = s/r;
  GV(3,3) = c/r;

  M4 = M3*GV;

  Vt = GV';

  if (abs(M4(2,2))<1.0e-08)

    c = M4(1,1);
    s = M4(1,2);
    r = sqrt(c^2 + s^2);

    GV = eye(3,3);
    GV(1,1) = c/r;
    GV(1,2) = -s/r;
    GV(2,1) = s/r;
    GV(2,2) = c/r;

    M5 = M4*GV;

    Vt = GV'*Vt;

%   Switch columns 2 and 3 
    P  = zeros(3,3);
    P(1,1) = 1;
    P(3,2) = 1;
    P(2,3) = 1;

    M5 = M5*P;

    Vt = P'*Vt;

%   Switch rows 2 and 3    
    P  = zeros(3,3);
    P(1,1) = 1;
    P(3,2) = 1;
    P(2,3) = 1;

    M5 = P*M5;

    Usvd  = Usvd*P';

  else

    disp('oh oh...')
    pause
  end

  M5 

end 

[u,s,v] = svd(M);



