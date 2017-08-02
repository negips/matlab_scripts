function [mass nek_mass DXM1 DYM1 DXM1D DYM1D RXM1 RYM1 SXM1 SYM1 convx convy convall convxd convyd convalld convxd_new convyd_new convalld_new Cx_fld Cy_fld gradm1x gradm1y gradm1xd gradm1yd intpm1d wtsvecd nek_conv lpx lpy lpall nek_lp lpbc forc nodal2spec2d x_coeff y_coeff Dx Dy w2m1 xm1 ym1 JACM1 JACM1D xm1d ym1d] = MESem2D5(Nx,Ny,Nxd,Nyd,xc,yc,ifplot)

%addpath '../../';

%    Definition of basis functions:
%    Definition of nodal values

%    U02-----U12-----U22
%     |       |       |
%     |       |       |
%    U01-----U11-----U21
%     |       |       |
%     |       |       |
%    U00-----U10-----U20

% For full matrix:
% Vector: {U00,
%          U10,
%          U20,
%          U01,
%          U11,
%          U21,
%          U02,
%          U12,
%          U22}

tic

epstol=1e-20;

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

w2m1 = zeros(Nx+1,Ny+1);
for i=0:Nx
for j=0:Ny
     w2m1(i+1,j+1) = wx(i+1)*wy(j+1);
end
end

N = Nx;                  % temporary

% Build polynomial coefficients for x 
A1 = [];
for i = 0:Nx
     A1= [A1 x.^i];
end
x_coeff = zeros(Nx+1);

for i = 0:Nx
     b = zeros(Nx+1,1);
     b(i+1) = 1;

     solns = A1\b;
     x_coeff(:,i+1) = solns;
end
%---------------------------------------- 
% Build polynomial coefficients for y 
A1 = [];
for i = 0:Ny
     A1= [A1 y.^i];
end
y_coeff = zeros(Ny+1);

for i = 0:Ny
     b = zeros(Ny+1,1);
     b(i+1) = 1;

     solns = A1\b;
     y_coeff(:,i+1) = solns;
end
%---------------------------------------- 
% Build derivative matrix (d/dx)
Dx = [];
for i = 0:Nx
     j = (i-1);
     t1 = i*x_coeff(i+1,:);
     nans = isnan(t1);
     if max(nans)
          ind = find(nans);
          t1(ind) = 0;
     end
     Dx= [Dx; t1];
end
%---------------------------------------- 
% Build derivative matrix (d/dy)
Dy = [];
for i = 0:Ny
     j = (i-1);
     t1 = i*y_coeff(i+1,:);
     nans = isnan(t1);
     if max(nans)
          ind = find(nans);
          t1(ind) = 0;
     end
     Dy= [Dy; t1];
end
%---------------------------------------- 

%% Plot Basis functions (in x) 
if ifplot
     xtemp = transpose(linspace(-1,1,200));
     l2 = length(xtemp);
     basis = zeros(l2,Nx+1);
     for i = 0:Nx
     val = zeros(l2,1);
     for j = 0:Nx
          val = val + x_coeff(j+1,i+1)*xtemp.^j;
     end
     basis(:,i+1) = val;
     end

     h1 = figure;
     hold on
     plot(xtemp,basis)
end

% Derivative function
if ifplot
     deriv = zeros(l2,Nx+1);
     for i = 0:Nx
     val = zeros(l2,1);
     for j = 0:Nx
          val1 = Dx(j+1,i+1)*(xtemp.^(j-1));
          if (j-1)<0
               val1 = zeros(l2,1);
          end
          val = val + val1;
     end
     deriv(:,i+1) = val;
     end

     h2 = figure;
     hold on
     plot(xtemp,deriv)
end
%---------------------------------------- 
% integral vector
dx_integral = zeros(2*N+1,1);
for j = 0:2*N
     xorder = j;
     dx_integral(j+1) = (1^(xorder+1) - (-1)^(xorder+1))/(xorder+1);
end
%---------------------------------------- 

disp('Calculating Jacobian')
%% Derivative operator d/dx
DXM1 = zeros(Nx+1,Nx+1);
d0 = Nx*(Nx+1)/4;
for ii=0:Nx
for jj=0:Nx
     DXM1(ii+1,jj+1) = 0;
     if ii~=jj
          lni = legendrePoly(Nx,x(ii+1));
          lnj = legendrePoly(Nx,x(jj+1));
          DXM1(ii+1,jj+1) = lni(end)/lnj(end)/(x(ii+1)-x(jj+1));
     end
     if ii==jj && ii==0
          DXM1(ii+1,jj+1) = -d0;
     end
     if ii==jj && ii==Nx
          DXM1(ii+1,jj+1) = d0;
     end
end
end
%---------------------------------------- 

%% Derivative operator d/dy
DYM1 = zeros(Ny+1,Ny+1);
d0 = Ny*(Ny+1)/4;
for ii=0:Ny
for jj=0:Ny
     DYM1(ii+1,jj+1) = 0;
     if ii~=jj
          lni = legendrePoly(Ny,y(ii+1));
          lnj = legendrePoly(Ny,y(jj+1));
          DYM1(ii+1,jj+1) = lni(end)/lnj(end)/(y(ii+1)-y(jj+1));
     end
     if ii==jj && ii==0
          DYM1(ii+1,jj+1) = -d0;
     end
     if ii==jj && ii==Ny
          DYM1(ii+1,jj+1) = d0;
     end
end
end
DYTM1 = transpose(DYM1);
%---------------------------------------- 

[xm1 ym1] = getlgll(Nx,Ny,xc,yc);

XRM1 = zeros(Nx+1,Ny+1);
YRM1 = zeros(Nx+1,Ny+1);
for ii=0:Ny
     XRM1(:,ii+1) = DXM1*xm1(:,ii+1);
     YRM1(:,ii+1) = DXM1*ym1(:,ii+1);
end

XSM1 = zeros(Nx+1,Ny+1);
YSM1 = zeros(Nx+1,Ny+1);
for ii=0:Nx
     XSM1(ii+1,:) = xm1(ii+1,:)*DYTM1;
     YSM1(ii+1,:) = ym1(ii+1,:)*DYTM1;
end

JACM1 = XRM1.*YSM1 - XSM1.*YRM1;

%---------------------------------------- 
% All these factors need to be divided by the jacobian (point wise multiplication).
% It get multiplied by the jacobian during integral.
% Hence it is skipped right now.
RXM1 = YSM1;
RYM1 = -XSM1;

SXM1 = -YRM1;
SYM1 = XRM1;

%---------------------------------------- 

% Mass matrix (in front of d/dt)
% This is a big but sparse matrix
disp('Calculating Mass matrix')
mass = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

for n=0:Ny
  for m=0:Nx
    ii=n*(Nx+1) + m+1;
    mass(ii,ii) = w2m1(m+1,n+1);
  end    
end    

%for n = 0:Ny
%     Ln = y_coeff(:,n+1);               % Test function (vx)
%for m = 0:Nx      
%     Lm = x_coeff(:,m+1);               % Test function (vx)
%     for j = 0:Ny
%          Lj = y_coeff(:,j+1);          % Trial function (uy)
%     for i = 0:Nx
%          Li = x_coeff(:,i+1);               % Trial Function (ux)
%
%%%        Integration
%          integral = 0;
%          for l =0:Ny                        % Weight along x
%               yl = y(l+1);
%          for k =0:Nx                        % Weight along y
%               xk = x(k+1);
%
%          %    Trial function/derivative values.
%               ifderiv =0;
%               Lm_xk = FuncEval(Lm,xk,ifderiv);
%               ifderiv =0;
%               Ln_yl = FuncEval(Ln,yl,ifderiv);
%
%          %    Test function/derivative values.
%               ifderiv =0;
%               Li_xk = FuncEval(Li,xk,ifderiv);
%               ifderiv =0;
%               Lj_yl = FuncEval(Lj,yl,ifderiv);
%
%               integral1 = JACM1(k+1,l+1)*w2m1(k+1,l+1)*Lm_xk*Ln_yl*(Li_xk*Lj_yl);
%               integral = integral + integral1; 
%          end
%          end
%%%        End of integration
%
%          posx = n*(Nx+1) + m + 1;
%          posy = j*(Nx+1) + i + 1;
%%         Since we build this very generally...
%%         Numerical errors occur
%          if abs(integral)<1e-12    
%               integral==0;
%          end
%          mass(posx,posy) = integral;
%     end
%     end
%
%end
%end

nek_mass = zeros(Nx+1,Ny+1);           % Only the diagonal terms are stored (as column vectors)
%
modeno=0;
for jj = 0:Ny     
     for ii = 0:Nx
          modeno=modeno+1;
          nek_mass(ii+1,jj+1) = mass(modeno,modeno);
     end
end


% matrix for vdu/dx.           
% This is a full matrix.
% Nek doesn't really build these matrices
% Making this in a very general way...
% disp('Calculating linear convective matrix')
 
convx = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

convy = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

convall = convx + convy;
%-------------------------------------------------- 
% Convective operator as NEK (might have) built.
nek_conv = zeros(Nx+1,Nx+1,Ny+1);

%--------------------------------------------------
 
%% Laplacian term                  % Integration by parts:  dudx.dvdx
disp('Calculating laplacian')
lpx = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

%for n = 0:Ny
%     Ln = y_coeff(:,n+1);               % Test function (vx)
%     dLn = Dy(:,n+1);
%for m = 0:Nx      
%     dLm = Dx(:,m+1);                        
%     Lm = x_coeff(:,m+1);                    % Test function (vx)
%     for j = 0:Ny
%          Lj = y_coeff(:,j+1);               % Trial function (uy)
%          dLj = Dy(:,j+1);
%     for i = 0:Nx
%          Li = x_coeff(:,i+1);               % Trial Function (ux)
%          dLi = Dx(:,i+1);
%
%%%        Integration
%          integral = 0;
%          for l =0:Ny                        % Weight along y
%               wty = wy(l+1);
%               yl = y(l+1);
%          for k =0:Nx                        % Weight along x
%               wtx = wx(k+1);
%               xk = x(k+1);
%
%          %    Trial function/derivative values.
%               ifderiv =0;
%               Lm_rk = FuncEval(Lm,xk,ifderiv);
%               ifderiv =0;
%               Ln_sl = FuncEval(Ln,yl,ifderiv);
%
%               ifderiv =1;
%               dLm_rk = FuncEval(dLm,xk,ifderiv);
%               ifderiv =1;
%               dLn_sl = FuncEval(dLn,yl,ifderiv);
%
%
%          %    Test function/derivative values.
%               ifderiv =0;
%               Li_rk = FuncEval(Li,xk,ifderiv);
%               ifderiv =0;
%               Lj_sl = FuncEval(Lj,yl,ifderiv);
%
%               ifderiv =1;
%               dLi_rk = FuncEval(dLi,xk,ifderiv);
%               ifderiv =1;
%               dLj_sl = FuncEval(dLj,yl,ifderiv);
%
%               integral1 = w2m1(k+1,l+1)*(dLm_rk*RXM1(k+1,l+1)*Ln_sl + ...
%                         dLn_sl*SXM1(k+1,l+1)*Lm_rk)* ...
%                         (dLi_rk*RXM1(k+1,l+1)*Lj_sl + dLj_sl*SXM1(k+1,l+1)*Li_rk);
%               integral = integral + integral1; 
%          end
%          end
%%%        End of integration
%
%          posx = n*(Nx+1) + m + 1;
%          posy = j*(Nx+1) + i + 1;
%%         Since we build this very generally...
%%         Numerical errors occur
%          if abs(integral)<1e-12
%               integral==0;
%          end
%          lpx(posx,posy) = integral;
%     end
%     end
%
%end
%end
%lpx = -lpx;
%
%
%% Laplacian term                  % Integration by parts:  dudx.dvdx
lpy = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
%
%for n = 0:Ny
%     Ln = y_coeff(:,n+1);                         % Test function (vx)
%     dLn = Dy(:,n+1);
%for m = 0:Nx      
%     Lm = x_coeff(:,m+1);                         % Test function (vx)
%     dLm = Dx(:,m+1);
%     for j = 0:Ny
%          Lj = y_coeff(:,j+1);                         % Trial function (uy)
%          dLj = Dy(:,j+1);
%     for i = 0:Nx
%          Li = x_coeff(:,i+1);                    % Trial Function (ux)
%          dLi = Dx(:,i+1);
%
%%%        Integration
%          integral = 0;
%          for l =0:Ny                             % Weight along y
%               wty = wy(l+1);
%               yl = y(l+1);
%          for k =0:Nx                             % Weight along x
%               wtx = wx(k+1);
%               xk = x(k+1);
%
%          %    Trial function/derivative values.
%               ifderiv =0;
%               Lm_rk = FuncEval(Lm,xk,ifderiv);
%               ifderiv =0;
%               Ln_sl = FuncEval(Ln,yl,ifderiv);
%
%               ifderiv =1;
%               dLm_rk = FuncEval(dLm,xk,ifderiv);
%               ifderiv =1;
%               dLn_sl = FuncEval(dLn,yl,ifderiv);
%
%
%          %    Test function/derivative values.
%               ifderiv =0;
%               Li_rk = FuncEval(Li,xk,ifderiv);
%               ifderiv =0;
%               Lj_sl = FuncEval(Lj,yl,ifderiv);
%
%               ifderiv =1;
%               dLi_rk = FuncEval(dLi,xk,ifderiv);
%               ifderiv =1;
%               dLj_sl = FuncEval(dLj,yl,ifderiv);
%
%               integral1 = w2m1(k+1,l+1)*(dLm_rk*RYM1(k+1,l+1)*Ln_sl + ...
%                         dLn_sl*SYM1(k+1,l+1)*Lm_rk)* ...
%                         (dLi_rk*RYM1(k+1,l+1)*Lj_sl + dLj_sl*SYM1(k+1,l+1)*Li_rk);
%               integral = integral + integral1;
%          end
%          end
%%%        End of integration
%
%          posx = n*(Nx+1) + m + 1;
%          posy = j*(Nx+1) + i + 1;
%%         Since we build this very generally...
%%         Numerical errors occur
%          if abs(integral)<1e-10
%               integral==0;
%          end
%          lpy(posx,posy) = integral;
%     end
%     end
%
%end
%end
%lpy = -lpy;
%
lpall = lpx + lpy;

%--------------------------------------------------
%    Laplacian operator as Nek builds it. 
nek_lp = zeros(Nx+1,Nx+1,Ny+1);
for kk=0:Ny
     ii_st = kk*(Nx+1)+1;
     ii_en = (kk+1)*(Nx+1);

     jj_st = kk*(Nx+1)+1;
     jj_en = (kk+1)*(Nx+1);
     nek_lp(:,:,kk+1) = lpall(ii_st:ii_en,jj_st:jj_en);
end
%-------------------------------------------------- 
lpbc = zeros(N+1);             % Boundary terms for laplacian operator. vdudx(1) - vdudx(-1)

%% Dealiasing
%-------------------------------------------------- 
disp('Calculating dealiased linear convective matrix')

if Nxd>1
     [xd wxd p]= lglnodes(Nxd);
     xd =xd(end:-1:1);
     wxd =wxd(end:-1:1);
else
     xd=[0];
     wxd=[1];
end

if Nyd>1
     [yd wyd p]= lglnodes(Nyd);
     yd =yd(end:-1:1);
     wyd = wyd(end:-1:1);
else
     yd=[0];
     wyd=[1];
end


% Build polynomial coefficients for xd 
A1 = [];
for i = 0:Nxd
     A1= [A1 xd.^i];
end
xd_coeff = zeros(Nxd+1);

for i = 0:Nxd
     b = zeros(Nxd+1,1);
     b(i+1) = 1;

     solns = A1\b;
     xd_coeff(:,i+1) = solns;
end
%---------------------------------------- 
% Build polynomial coefficients for yd 
A1 = [];
for i = 0:Nyd
     A1= [A1 yd.^i];
end
yd_coeff = zeros(Nyd+1);

for i = 0:Nyd
     b = zeros(Nyd+1,1);
     b(i+1) = 1;

     solns = A1\b;
     yd_coeff(:,i+1) = solns;
end

%---------------------------------------- 
% Build derivative matrix (d/dx)
Dxd = [];
for i = 0:Nxd
     j = (i-1);
     t1 = i*xd_coeff(i+1,:);
     nans = isnan(t1);
     if max(nans)
          ind = find(nans);
          t1(ind) = 0;
     end
     Dxd= [Dxd; t1];
end
%---------------------------------------- 
% Build derivative matrix (d/dy)
Dyd = [];
for i = 0:Nyd
     j = (i-1);
     t1 = i*yd_coeff(i+1,:);
     nans = isnan(t1);
     if max(nans)
          ind = find(nans);
          t1(ind) = 0;
     end
     Dyd= [Dyd; t1];
end
%---------------------------------------- 


%% Derivative operator d/dx
DXM1D = zeros(Nxd+1,Nx+1);
for ii=0:Nxd
     xi = xd(ii+1);
for jj=0:Nx
     dLj = Dx(:,jj+1);
     ifderiv=1;
     dLj_xi = FuncEval(dLj,xi,ifderiv);

     DXM1D(ii+1,jj+1) = dLj_xi;
end
end

%---------------------------------------- 

%% Derivative operator d/dy
DYM1D = zeros(Nyd+1,Ny+1);
for ii=0:Nyd
     yi = yd(ii+1);
for jj=0:Ny
     dLj = Dy(:,jj+1);
     ifderiv=1;
     dLj_yi = FuncEval(dLj,yi,ifderiv);

     DYM1D(ii+1,jj+1) = dLj_yi;
end
end

DYTM1D = transpose(DYM1D);
%----------------------------------------
[xm1d ym1d] = getlgll(Nxd,Nyd,xc,yc);

XRM1D = zeros(Nxd+1,Nyd+1);
YRM1D = zeros(Nxd+1,Nyd+1);
for ii=0:Nyd

     xtmp = getlgll1D(Nx,xm1d(1,ii+1),xm1d(end,ii+1));
     ytmp = getlgll1D(Nx,ym1d(1,ii+1),ym1d(end,ii+1));

     XRM1D(:,ii+1) = DXM1D*xtmp;
     YRM1D(:,ii+1) = DXM1D*ytmp;
end

XSM1D = zeros(Nxd+1,Nyd+1);
YSM1D = zeros(Nxd+1,Nyd+1);
for ii=0:Nxd

     xtmp = getlgll1D(Ny,xm1d(ii+1,1),xm1d(ii+1,end));
     ytmp = getlgll1D(Ny,ym1d(ii+1,1),ym1d(ii+1,end));

     XSM1D(ii+1,:) = transpose(xtmp)*DYTM1D;
     YSM1D(ii+1,:) = transpose(ytmp)*DYTM1D;
end

JACM1D = XRM1D.*YSM1D - XSM1D.*YRM1D;

%---------------------------------------- 
% All these factors need to be divided by the jacobian (point wise multiplication).
% It get multiplied by the jacobian during integral.
% Hence it is skipped right now.
RXM1D = YSM1D;
RYM1D = -XSM1D;

SXM1D = -YRM1D;
SYM1D = XRM1D;

% matrix for vdu/dx.           
% This is a full matrix.
% Nek doesn't really build these matrices
% Making this in a very general way...
[Cx_fld Cy_fld] = GetConvectingfield(xm1,ym1);       % On the GLL grid
[Cxd_fld Cyd_fld] = GetConvectingfield(xm1d,ym1d);       % On the dealiased grid


convxd = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
%
%for n = 0:Ny
%     Ln = y_coeff(:,n+1);               % Test function (vx)
%for m = 0:Nx      
%     Lm = x_coeff(:,m+1);               % Test function (vx)
%     for j = 0:Ny
%          Lj = y_coeff(:,j+1);          % Trial function (uy)
%          dLj = Dy(:,j+1);
%     for i = 0:Nx
%          Li = x_coeff(:,i+1);               % Trial Function (ux)
%          dLi = Dx(:,i+1);
%
%%%        Integration
%          integral = 0;
%          for k =0:Nxd                        % Weight along x
%               wtx = wxd(k+1);
%               xk = xd(k+1);
%          for l =0:Nyd                        % Weight along y
%               wty = wyd(l+1);
%               yl = yd(l+1);
%
%          %    Trial function/derivative values.
%               ifderiv =0;
%               Lm_xk = FuncEval(Lm,xk,ifderiv);
%               ifderiv =0;
%               Ln_yl = FuncEval(Ln,yl,ifderiv);
%
%          %    Test function/derivative values.
%               ifderiv =0;
%               Li_xk = FuncEval(Li,xk,ifderiv);
%               ifderiv =0;
%               Lj_yl = FuncEval(Lj,yl,ifderiv);
%
%               ifderiv =1;
%               dLi_xk = FuncEval(dLi,xk,ifderiv);
%               ifderiv =1;
%               dLj_yl = FuncEval(dLj,yl,ifderiv);
%
%%              We don't multiply with jacobian here since RXM1 etc have (JACM1)^(-1) factor.
%%              Which is omitted in those terms
%               integral1 = wtx*wty*Lm_xk*Ln_yl*(dLi_xk*RXM1D(k+1,l+1)*Lj_yl*Cxd_fld(k+1,l+1) + Li_xk*SXM1D(k+1,l+1)*dLj_yl*Cyd_fld(k+1,l+1));
%               integral = integral + integral1; 
%          end
%          end
%%%        End of integration
%
%          posx = n*(Nx+1) + m + 1;
%          posy = j*(Nx+1) + i + 1;
%%         Since we build this very generally...
%%         Numerical errors occur
%          if abs(integral)<1e-12
%               integral==0;
%          end
%          convxd(posx,posy) = integral;
%     end
%     end
%
%end
%end
%
convyd = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
%
%for n = 0:Ny
%     Ln = y_coeff(:,n+1);               % Test function (vx)
%for m = 0:Nx      
%     Lm = x_coeff(:,m+1);               % Test function (vx)
%     for j = 0:Ny
%          Lj = y_coeff(:,j+1);          % Trial function (uy)
%          dLj = Dy(:,j+1);
%     for i = 0:Nx
%          Li = x_coeff(:,i+1);               % Trial Function (ux)
%          dLi = Dx(:,i+1);
%
%%%        Integration
%          integral = 0;
%          for k =0:Nxd                        % Weight along x
%               wtx = wxd(k+1);
%               xk = xd(k+1);
%          for l =0:Nyd                        % Weight along y
%               wty = wyd(l+1);
%               yl = yd(l+1);
%
%          %    Trial function/derivative values.
%               ifderiv =0;
%               Lm_xk = FuncEval(Lm,xk,ifderiv);
%               ifderiv =0;
%               Ln_yl = FuncEval(Ln,yl,ifderiv);
%
%          %    Test function/derivative values.
%               ifderiv =0;
%               Li_xk = FuncEval(Li,xk,ifderiv);
%               ifderiv =0;
%               Lj_yl = FuncEval(Lj,yl,ifderiv);
%
%               ifderiv =1;
%               dLi_xk = FuncEval(dLi,xk,ifderiv);
%               ifderiv =1;
%               dLj_yl = FuncEval(dLj,yl,ifderiv);
%
%%              We don't multiply with jacobian here since RXM1 etc have (JACM1)^(-1) factor.
%%              Which is omitted in those terms
%               integral1 = wtx*wty*Lm_xk*Ln_yl*(dLi_xk*RYM1D(k+1,l+1)*Lj_yl + Li_xk*SYM1D(k+1,l+1)*dLj_yl);
%               integral = integral + integral1; 
%
%          end
%          end
%%%        End of integration
%
%          posx = n*(Nx+1) + m + 1;
%          posy = j*(Nx+1) + i + 1;
%%         Since we build this very generally...
%%         Numerical errors occur
%          if abs(integral)<1e-12
%               integral==0;
%          end
%          convyd(posx,posy) = integral;
%     end
%     end
%
%end
%end
%
convalld = convxd + convyd;

% Just derivative operator.
% No integration
% du/dx (on the dealiased grid).
gradm1x = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

%%        De-aliasing points 
for l =0:Ny                        % points along s
     wts = wy(l+1);
     sl = y(l+1);

for k =0:Nx                        % points along r
     wtr = wx(k+1);
     rk = x(k+1);

     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (us)
          dLj = Dy(:,j+1);
     for i = 0:Nx
          Li = x_coeff(:,i+1);          % Trial Function (ur)
          dLi = Dx(:,i+1);

%         Test function/derivative values.
          ifderiv =0;
          Li_rk = FuncEval(Li,rk,ifderiv);
          ifderiv =0;
          Lj_sl = FuncEval(Lj,sl,ifderiv);

          ifderiv =1;
          dLi_rk = FuncEval(dLi,rk,ifderiv);
          ifderiv =1;
          dLj_sl = FuncEval(dLj,sl,ifderiv);

          deriv_val = (dLi_rk*RXM1(k+1,l+1)*Lj_sl + Li_rk*dLj_sl*SXM1(k+1,l+1))/JACM1(k+1,l+1);

          posx = l*(Nx+1) + k + 1;
          posy = j*(Nx+1) + i + 1;
%         Since we build this very generally...
%         Numerical errors occur (O(1e-14))
          if (abs(deriv_val)>epstol)  
            gradm1x(posx,posy) = deriv_val;
          end  

     end
     end

end
end

% Just derivative operator.
% No integration
% du/dy (on the dealiased grid).

gradm1y = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

%%        De-aliasing points 
for l =0:Ny                        % points along s
     wts = wy(l+1);
     sl = y(l+1);

for k =0:Nx                        % points along r
     wtr = wx(k+1);
     rk = x(k+1);

     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (us)
          dLj = Dy(:,j+1);
     for i = 0:Nx
          Li = x_coeff(:,i+1);          % Trial Function (ur)
          dLi = Dx(:,i+1);

%         Test function/derivative values.
          ifderiv =0;
          Li_rk = FuncEval(Li,rk,ifderiv);
          ifderiv =0;
          Lj_sl = FuncEval(Lj,sl,ifderiv);

          ifderiv =1;
          dLi_rk = FuncEval(dLi,rk,ifderiv);
          ifderiv =1;
          dLj_sl = FuncEval(dLj,sl,ifderiv);

          deriv_val = (dLi_rk*RYM1(k+1,l+1)*Lj_sl + Li_rk*dLj_sl*SYM1(k+1,l+1))/JACM1(k+1,l+1);

          posx = l*(Nx+1) + k + 1;
          posy = j*(Nx+1) + i + 1;
%         Since we build this very generally...
%         Numerical errors occur (O(1e-14))
          if (abs(deriv_val)>epstol)  
            gradm1y(posx,posy) = deriv_val;
          end  

     end
     end

end
end


% Just derivative operator.
% No integration
% du/dx (on the dealiased grid).

disp('Calculating derivative operators')
gradm1xd = zeros((Nxd+1)*(Nyd+1),(Nx+1)*(Ny+1));

%%        De-aliasing points 
for l =0:Nyd                        % points along s
     wts = wyd(l+1);
     sl = yd(l+1);

for k =0:Nxd                        % points along r
     wtr = wxd(k+1);
     rk = xd(k+1);

     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (us)
          dLj = Dy(:,j+1);
     for i = 0:Nx
          Li = x_coeff(:,i+1);          % Trial Function (ur)
          dLi = Dx(:,i+1);

%         Test function/derivative values.
          ifderiv =0;
          Li_rk = FuncEval(Li,rk,ifderiv);
          ifderiv =0;
          Lj_sl = FuncEval(Lj,sl,ifderiv);

          ifderiv =1;
          dLi_rk = FuncEval(dLi,rk,ifderiv);
          ifderiv =1;
          dLj_sl = FuncEval(dLj,sl,ifderiv);

          deriv_val = (dLi_rk*RXM1D(k+1,l+1)*Lj_sl + Li_rk*dLj_sl*SXM1D(k+1,l+1))/JACM1D(k+1,l+1);

          posx = l*(Nxd+1) + k + 1;
          posy = j*(Nx+1) + i + 1;
%         Since we build this very generally...
%         Numerical errors occur (O(1e-14))
          if (abs(deriv_val)>epstol)  
            gradm1xd(posx,posy) = deriv_val;
          end  


     end
     end

end
end

% Just derivative operator.
% No integration
% du/dy (on the dealiased grid).

gradm1yd = zeros((Nxd+1)*(Nyd+1),(Nx+1)*(Ny+1));

%%        De-aliasing points 
for l =0:Nyd                        % points along s
     wts = wyd(l+1);
     sl = yd(l+1);

for k =0:Nxd                        % points along r
     wtr = wxd(k+1);
     rk = xd(k+1);

     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (us)
          dLj = Dy(:,j+1);
     for i = 0:Nx
          Li = x_coeff(:,i+1);          % Trial Function (ur)
          dLi = Dx(:,i+1);

%         Test function/derivative values.
          ifderiv =0;
          Li_rk = FuncEval(Li,rk,ifderiv);
          ifderiv =0;
          Lj_sl = FuncEval(Lj,sl,ifderiv);

          ifderiv =1;
          dLi_rk = FuncEval(dLi,rk,ifderiv);
          ifderiv =1;
          dLj_sl = FuncEval(dLj,sl,ifderiv);

          deriv_val = (dLi_rk*RYM1D(k+1,l+1)*Lj_sl + Li_rk*dLj_sl*SYM1D(k+1,l+1))/JACM1D(k+1,l+1);

          posx = l*(Nxd+1) + k + 1;
          posy = j*(Nx+1) + i + 1;
%         Since we build this very generally...
%         Numerical errors occur (O(1e-14))
          if (abs(deriv_val)>epstol)  
            gradm1yd(posx,posy) = deriv_val;
          end  
     end
     end

end
end


intpm1d = zeros((Nxd+1)*(Nyd+1),(Nx+1)*(Ny+1));
wtsvecd = zeros((Nxd+1)*(Nyd+1),1);       % Weight vector for dealiased grid.
%%        De-aliasing points 
for l =0:Nyd                        % points along y
     wts = wyd(l+1);
     sl = yd(l+1);

for k =0:Nxd                        % points along x
     wtr = wxd(k+1);
     rk = xd(k+1);

     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (uy)
          dLj = Dy(:,j+1);
     for i = 0:Nx
          Li = x_coeff(:,i+1);          % Trial Function (ux)
          dLi = Dx(:,i+1);

%         Test function/derivative values.
          ifderiv =0;
          Li_rk = FuncEval(Li,rk,ifderiv);
          ifderiv =0;
          Lj_sl = FuncEval(Lj,sl,ifderiv);

%         No derivates here ...

          intp_val = Li_rk*Lj_sl;

          posx = l*(Nxd+1) + k + 1;
          posy = j*(Nx+1) + i + 1;
%         Since we build this very generally...
%         Numerical errors occur (O(1e-14))
          intpm1d(posx,posy) = intp_val;

     end
     end

     posx = l*(Nxd+1) + k + 1;
     wtsvecd(posx,1) = wtr*wts;

end
end


%% Build forcing

%% Build basis change matrix
%% Taken directly from nek
%% In the x-direction
nx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:nx
  z=x(j);
  Lj = legendrePoly(n,z);
  kj = kj+1;
  pht(kj) = Lj(1);
  kj = kj+1;
  pht(kj) = Lj(2);
  for k = 3:nx
    kj = kj+1;
    pht(kj) = Lj(k);% - Lj(k-2);
  end
end

pht = reshape(pht,nx,nx);
spec2nodalx = transpose(pht);
nodal2specx = inv(spec2nodalx);

% Filter function. Modes and amplitudes
alpha_x=0.25;
Gx = eye(Nx+1);
Gx(Nx+1,Nx+1)=1-alpha_x;
%Gx(Nx,Nx)=1-alpha_x.^2;

% Build matrix which applies filter on Boyd transformed basis.
% And then transforms back to physical space.
phi = pht;

boyd = 0;
if boyd
  for i = 3:nx
    for j = 1:nx
      phi(i,j) = pht(i,j) - pht(i-2,j);
    end
  end
end
boydstonodal_x = transpose(phi);
nodaltoboyds_x = inv(boydstonodal_x);

fil_mat_x = boydstonodal_x*Gx*nodaltoboyds_x;


%% In the y-direction
ny = Ny+1;
kj = 0;
n = ny-1;
for j = 1:ny
  z=x(j);
  Lj = legendrePoly(n,z);
  kj = kj+1;
  pht(kj) = Lj(1);
  kj = kj+1;
  pht(kj) = Lj(2);
  for k = 3:ny
    kj = kj+1;
    pht(kj) = Lj(k);% - Lj(k-2);
  end
end

pht = reshape(pht,ny,ny);
spec2nodaly = transpose(pht);
nodal2specy = inv(spec2nodaly);

phi = pht;

boyd = 0;
if boyd
  for i = 3:ny
    for j = 1:ny
      phi(i,j) = pht(i,j) - pht(i-2,j);
    end
  end
end
boydstonodal_y = transpose(phi);
nodaltoboyds_y = inv(boydstonodal_y);


% Filter function. Modes and amplitudes
alpha_y=0.25;
Gy = eye(Ny+1);
Gy(Ny+1,Ny+1)=1-alpha_y;
%Gy(Ny,Ny)=1-alpha_y.^2;

% Build matrix which applies filter on Boyd transformed basis.
% And then transforms back to physical space.
fil_mat_y = boydstonodal_y*Gy*nodaltoboyds_y;
fil_mat2d = kron(fil_mat_x,fil_mat_y);

spec2nodal2d = kron(spec2nodalx,spec2nodaly);
nodal2spec2d = kron(nodal2specx,nodal2specy);


forc = mass*(eye((Nx+1)*(Ny+1))-fil_mat2d);


%% Build reverse interpolation
% Dealiasing grid to GLL grid

%% Build dealiasing matrix
%% In the x-direction
pht = zeros(Nx+1,Nxd+1);
nx = Nxd+1;
lx = Nx+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=x(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
  if (Nxd>Nx)    
    pht(j,lx+1:nx) = 0;               % Truncate space to order N
  end 
end

dealiased_intpx = pht;

%% In the y-direction
pht = zeros(Ny+1,Nyd+1);
ny = Nyd+1;
ly = Ny+1;
kj = 0;
n = ny-1;
for j = 1:ly
  z=y(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
  if (Nyd>Ny)
    pht(j,ly+1:ny) = 0;               % Truncate space to order N
  end    
end

dealiased_intpy = pht;

dealias_spec2gll = kron(dealiased_intpx,dealiased_intpy);

%% Build nodal to spectral transform in dealiasing grid
pht = zeros(Nxd+1,Nxd+1);
nx = Nxd+1;
lx = Nxd+1;
kj = 0;
n = nx-1;
for j = 1:lx
  z=xd(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end

dealiased_nodal2specx = inv(pht);

%% In the y-direction
pht = zeros(Nyd+1,Nyd+1);
ny = Nyd+1;
ly = Nyd+1;
kj = 0;
n = ny-1;
for j = 1:ly
  z=yd(j);
  Lj = legendrePoly(n,z);
  pht(j,:) = transpose(Lj);
end

dealiased_nodal2specy = inv(pht);

dealiased_nodal2spec = kron(dealiased_nodal2specx,dealiased_nodal2specy);

dealias2gll = dealias_spec2gll*dealiased_nodal2spec;

%% Build fast dealiased-convection operator
%convxd_new = mass*dealias2gll*diag(Cxd_fld(:))*gradm1xd;
%convyd_new = mass*dealias2gll*diag(Cyd_fld(:))*gradm1yd;
%convalld_new = convxd_new + convyd_new;

convxd_new = mass*dealias2gll*diag(intpm1d*Cx_fld(:))*gradm1xd;
convyd_new = mass*dealias2gll*diag(intpm1d*Cy_fld(:))*gradm1yd;
convalld_new = convxd_new + convyd_new;



%dbstop in MESem2D6 at 1168

toc

return



