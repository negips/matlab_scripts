%  Build the big mass and stiffness matrices

   clear
   clc
%   close all

   spec_element_init_ms

   Umax=U+2*cd*cu;
%   Abs_limit = (Umax^2)/(4*abs(gamma)^2);
%   disp(['Absolute instability Limit, mu00=',num2str(Abs_limit,5) ])
%   if (mu00>Abs_limit)
%     disp(['Globally unstable, mu00=',num2str(mu00,5) ])
%   elseif (mu00<Abs_limit && mu00>0)
%     disp(['Convectively unstable, mu00=',num2str(mu00,5) ])
%   elseif (mu00<0)
%     disp(['Stable',num2str(mu00,5) ])
%   else
%     disp(['Marginaly Globaly Stable, mu00=',num2str(mu00,5) ])
%   end  

   gl_mass=zeros(dof,dof);
   gl_conv=zeros(dof,dof);
   gl_lapl=zeros(dof,dof);
   gl_src =zeros(dof,dof);

   for els=1:nels

      gl_pos_j1 = (els-1)*N + 1;
      gl_pos_i1 = (els-1)*N + 1;
  
      gl_pos_j2 = (els)*N + 1;
      gl_pos_i2 = (els)*N + 1;


      gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = nek_mass(:,:,els) + gl_mass(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
     
%     \nu*d/dx     
      gl_conv(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = nek_conv(:,:,els) + gl_conv(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);

%     \gamma*d^2/dx^2
      gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = nek_lp(:,:,els) + gl_lapl(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);

%     \mu
      gl_src(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2) = nek_mass(:,:,els).*nek_mu(:,:,els) + gl_src(gl_pos_j1:gl_pos_j2,gl_pos_i1:gl_pos_i2);
    
   end

   gl_stiff = -gl_conv - gl_lapl + gl_src;

%  Remove first boundary point    
   gl_stiff(1,:) = [];
   gl_stiff(:,1) = [];

%  Remove first boundary point    
   gl_mass(1,:) = [];
   gl_mass(:,1) = [];


%%  Remove last boundary point    
%   gl_stiff(end,:) = [];
%   gl_stiff(:,end) = [];
%
%%  Remove first boundary point    
%   gl_mass(end,:) = [];
%   gl_mass(:,end) = [];


    [V,D]=eig(gl_stiff,gl_mass);

    eg = diag(D);
    er = real(eg);
    ei = imag(eg);

    scatter(er,ei, '.r'); hold on
    grid on
