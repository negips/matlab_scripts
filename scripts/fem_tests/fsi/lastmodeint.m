function [mass] = lastmodeint(Nx,Ny,xc,yc)

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

% Mass matrix (in front of d/dt)
% This is a big but sparse matrix
%disp('Calculating Mass matrix')
mass = zeros((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));

for n = 0:Ny
     Ln = y_coeff(:,n+1);               % Test function (vx)
for m = 0:Nx      
     Lm = x_coeff(:,m+1);               % Test function (vx)
     for j = 0:Ny
          Lj = y_coeff(:,j+1);          % Trial function (uy)
     for i = 0:Nx
          Li = x_coeff(:,i+1);               % Trial Function (ux)

%%        Integration
          integral = 0;
          for l =0:Ny                        % Weight along x
               yl = y(l+1);
          for k =0:Nx                        % Weight along y
               xk = x(k+1);

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

               integral1 = JACM1(k+1,l+1)*w2m1(k+1,l+1)*Lm_xk*Ln_yl*(Li_xk*Lj_yl);
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
          mass(posx,posy) = integral;
     end
     end

end
end


return



