function [massd] = lastmodeintd(Nx,Ny,Nxd,Nyd,xc,yc)

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

%disp('Calculating Jacobian')
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

%% Dealiasing
%-------------------------------------------------- 
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

% dealiased matrix.
% This is a full matrix.
% Making this in a very general way...
massd = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

for n = 0:Ny
     Ln = y_coeff(:,n+1);               % Test function (vx)
for m = 0:Nx      
     Lm = x_coeff(:,m+1);               % Test function (vx)
     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (uy)
          dLj = Dy(:,j+1);
     for i = 0:Nx
          Li = x_coeff(:,i+1);               % Trial Function (ux)
          dLi = Dx(:,i+1);

%%        Integration
          integral = 0;
          for k =0:Nxd                        % Weight along x
               wtx = wxd(k+1);
               xk = xd(k+1);
          for l =0:Nyd                        % Weight along y
               wty = wyd(l+1);
               yl = yd(l+1);

          %    Trial function/derivative values.
               ifderiv =0;
               Lm_xk = FuncEval(Lm,xk,ifderiv);
               ifderiv =0;
               Ln_yl = FuncEval(Ln,yl,ifderiv);

          %    Test function/derivative values.
               ifderiv =0;
               Li_xk = FuncEval(Li,xk,ifderiv);
               ifderiv =0;
               Lj_yl = FuncEval(Lj,yl,ifderiv);

               ifderiv =1;
               dLi_xk = FuncEval(dLi,xk,ifderiv);
               ifderiv =1;
               dLj_yl = FuncEval(dLj,yl,ifderiv);

               integral1 = wtx*wty*Lm_xk*Ln_yl*(Li_xk*Lj_yl)*JACM1D(k+1,l+1);
               integral = integral + integral1; 
          end
          end
%%        End of integration

          posx = n*(Nx+1) + m + 1;
          posy = j*(Nx+1) + i + 1;
%         Since we build this very generally...
%         Numerical errors occur
          if abs(integral)<1e-12
               integral==0;
          end
          massd(posx,posy) = integral;
     end
     end

end
end

return



