%    Testing integration with mass matricies

function [int_fld] = dealias_int(Nx,Nxd1,Nxd2) 

%     Nx   = 6;
%     Nxd1 = 11;       % dealiasing points.
%     Nxd2 = 8;

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


%    dealiased product of 2 terms 
     vw_fine = zeros(Nxd1+1,1);

     for i=1:Nxd1+1
          xpt=xd1(i);

          Pn = legendrePoly(Nx,xpt);
          vw_fine(i) = Pn(end)*Pn(end);
     end

     [M x_coeff] = mass_mat(Nx);
%% Plot Basis functions (in x)
     ifplot=0; 
     if ifplot
%          xtemp = transpose(linspace(-1,1,500));
          xtemp = lglnodes(500);
          xtemp = xtemp(end:-1:1);
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
%-------------------- 

%    Dealiased integration
     int_fld = zeros(Nx+1,1); 
     for m=1:Nx+1
          Lm = x_coeff(:,m);

          integral = 0;
          for k=1:Nxd1+1

               xpt = xd1(k);
               wts = wxd1(k);

               ifderiv =0;
               Lm_xk = FuncEval(Lm,xpt,ifderiv);
               nodal_val = Lm_xk*vw_fine(k,1);
               integral = integral + wts*nodal_val;

          end
          int_fld(m) = integral;

     end

%     [spectonodal nodaltospec] = Leg2Nodal(Nxd1);
%     leg_fine = nodaltospec*vw_fine;
%
%     proj_coarse = proj_op(Nx,Nxd1);
%     nodal_coarse = proj_coarse*leg_fine;
%     integral = M*nodal_coarse;
     
end

%-------------------------------------------------- 

function A = proj_op(Nx,Nxd)

     x=lglnodes(Nx);
     x=x(end:-1:1);

     A = zeros(Nx+1,Nxd+1);
     for i = 1:Nx+1
          A(i,:) = transpose(legendrePoly(Nxd,x(i)));
     end
end

%--------------------------------------------------
 
function [M x_coeff] = mass_mat(Nx)

     [x wx p]=lglnodes(Nx);
     x=x(end:-1:1);
     wx = wx(end:-1:1);

%    Build polynomial coefficients for x 
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

     M = zeros(Nx+1,Nx+1);

     for m = 0:Nx      
          Lm = x_coeff(:,m+1);                    % Test function (vx)
          for i = 0:Nx
               Li = x_coeff(:,i+1);               % Trial Function (ux)

%%             Integration
               integral = 0;
               for k =0:Nx                        % Weight along x
                    xk = x(k+1);

               %    Trial function/derivative values.
                    ifderiv =0;
                    Lm_xk = FuncEval(Lm,xk,ifderiv);

               %    Test function/derivative values.
                    ifderiv =0;
                    Li_xk = FuncEval(Li,xk,ifderiv);

                    integral1 = wx(k+1)*Lm_xk*Li_xk;
                    integral = integral + integral1; 
               end
%%             End of integration

%              Since we build this very generally...
%              Numerical errors occur
               if abs(integral)<1e-12    
                    integral==0;
               end
               M(m+1,i+1) = integral;
          end
     end

%%

end

%-------------------------------------------------- 



