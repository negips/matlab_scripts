%    Testing integration with mass matricies

function [int_fld1 spec_l1 int_fld2 spec_l2] = paul_int(Nx,Nxd1,Nxd2) 

%     Nx   = 20;
%     Nxd1 = 21;       % dealiasing points.
%     Nxd2 = 8;

%    Calculating uvwdx 
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

     % dealiased product of 2 terms. 
     vw_fine = zeros(Nxd1+1,1);
     for i=1:Nxd1+1
          xpt=xd1(i);

          Pn = legendrePoly(Nx,xpt);
          vw_fine(i) = Pn(end)*Pn(end);
     end


     [spectonodal nodaltospec] = Leg2Nodal(Nxd1);
     leg_fine = nodaltospec*vw_fine;

     spec_l1 = leg_fine(1:Nx+1);
     leg_fine(Nx+2:end) = 0;

     proj_coarse = proj_op(Nx,Nxd1);
     vw_coarse = proj_coarse*leg_fine;

     [M x_coeff] = mass_mat(Nx);
     
     int_fld1 = M*vw_coarse;

     if Nxd2>0
          [xd wxd p]= lglnodes(Nxd2);
          xd2 =xd(end:-1:1);
          wxd2 =wxd(end:-1:1);
     end

     % dealiased product of 2 terms. 
     vw_fine = zeros(Nxd2+1,1);
     for i=1:Nxd2+1
          xpt=xd2(i);

          Pn = legendrePoly(Nx,xpt);
          vw_fine(i) = Pn(end)*Pn(end-1);
     end

     [spectonodal nodaltospec] = Leg2Nodal(Nxd2);
     leg_fine = nodaltospec*vw_fine;

     spec_l2 = leg_fine(1:Nx+1);
     leg_fine(Nx+2:end) = 0;

     proj_coarse = proj_op(Nx,Nxd2);
     vw_coarse = proj_coarse*leg_fine;

     [M x_coeff] = mass_mat(Nx);
     
     int_fld2 = M*vw_coarse;

    
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



