function c_op = convop(Nx,Ny,Nxd,Nyd,El,nelv,nbasis,ifnek)

Eltmp = [];

for el=1:nelv

     Eltmp(el).dudx = El(el).gradm1xd*El(el).un(:);
     Eltmp(el).dudy = El(el).gradm1yd*El(el).un(:);

     Eltmp(el).cxd = El(el).intpm1d*El(el).cx(:);
     Eltmp(el).cyd = El(el).intpm1d*El(el).cy(:);

end

if (ifnek)
     El = integration_nek(Eltmp,El,Nx,Ny,Nxd,Nyd,nelv,nbasis);
else
     El = integration_dealias(Eltmp,El,Nx,Ny,Nxd,Nyd,nelv,nbasis);
end

gl_rhs = DSSUMfld(El,Nx,Ny,nelv,nbasis);
c_op = gl_rhs;

%-------------------------------------------------- 


function Elout = integration_dealias(Eltmp,El,Nx,Ny,Nxd,Nyd,nelv,nbasis)

lbasis = (Nx+1)*(Ny+1);
rhs = zeros(lbasis,1);

for el=1:nelv

     for test = 1:lbasis
          
          test_fn = zeros(lbasis,1);
          test_fn(test) = 1;
          test_fnd = El(el).intpm1d*test_fn;

          tri_prod = (Eltmp(el).dudx.*Eltmp(el).cxd + Eltmp(el).dudy.*Eltmp(el).cyd).*test_fnd;

%         integration 
          integ = El(el).wtsvecd.*tri_prod.*El(el).JACM1D(:);
          rhs(test) = sum(integ);
     end

     El(el).local_fld = reshape(rhs,Nx+1,Ny+1);
end
Elout = El;
          
%-------------------------------------------------- 


function Elout = integration_nek(Eltmp,El,Nx,Ny,Nxd,Nyd,nelv,nbasis)
     Elout = El;

%-------------------------------------------------- 


function gl_fld = DSSUMfld(El,Nx,Ny,nelv,nbasis)

     gl_fld = zeros(nbasis,1);

     for el = 1:nelv
          for i=1:Nx+1
          for j=1:Ny+1
               pos = El(el).lglmap(i,j);
               gl_fld(pos) = gl_fld(pos) + El(el).local_fld(i,j);

          end
          end
     end
         
%-------------------------------------------------- 

