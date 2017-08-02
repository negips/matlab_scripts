function [dt dt2 dc dc2 VCdx lp cbc lpbc forc x_coeffs D x wts nodaltospec spectonodal] = MEFem_Cx(N,Nxd,xs,xe,ifplot)

jac = (xe-xs)/2;              % jacobian
dxbydxi = (xe-xs)/2;
dxibydx = 2/(xe-xs);

% GLL grid
[x wts p]= lglnodes(N);
x =x(end:-1:1);
wts =wts(end:-1:1);

% Dealiased grid
[xd wtsd pd]= lglnodes(Nxd);
xd =xd(end:-1:1);
wtsd =wtsd(end:-1:1);


%% build A
A1 = [];
for i = 0:N
  A1= [A1 x.^i];
end

x_coeffs = zeros(N+1);
for i = 0:N
  % build b
  b = zeros(N+1,1);
  b(i+1) = 1;

  solns = A1\b;
  x_coeffs(:,i+1) = solns;
end


%% Derivative
D = [];
for i = 0:N
  j = (i-1);
  t1 = i*x_coeffs(i+1,:);
  nans = isnan(t1);
  if max(nans)
    ind = find(nans);
    t1(ind) = 0;
  end
  D= [D; t1];
end


%x_coeffs
%D

%% testing
xtemp = transpose(linspace(-1,1,200));
l2 = length(xtemp);
basis = zeros(l2,N+1);
for i = 0:N
  val = zeros(l2,1);
  for j = 0:N
    val = val + x_coeffs(j+1,i+1)*xtemp.^j;
  end
  basis(:,i+1) = val;
end

if ifplot
  h1 = figure;
  hold on
  plot(xtemp,basis)
end

% Derivative values on GLL grid
deriv = zeros(l2,N+1);
for i = 0:N
  val = zeros(l2,1);
  for j = 0:N
    val1 = D(j+1,i+1)*(xtemp.^(j-1));
    if (j-1)<0
       val1 = zeros(l2,1);
    end
    val = val + val1;
  end
  deriv(:,i+1) = val;
end

% Derivative values on dealiased grid
deriv_d = zeros(l2,N+1);
xtmp_deriv = zeros(N+1,l2);
for i = 0:N
  if (i==0)
    xtmp_deriv(i+1,:) = zeros(1,l2);
  else
    xtmp_deriv(i+1,:)=transpose(xtemp.^(i-1));
  end
end
deriv_2 = transpose(transpose(D)*xtmp_deriv);


if ifplot
  h2 = figure;
  hold on
  plot(xtemp,deriv); hold on
  plot(xtemp,deriv_2, '--')
end
%% End of testing


% integral vector
dx_integral = zeros(2*N+1,1);
for j = 0:2*N
  xorder = j;
  dx_integral(j+1) = (1^(xorder+1) - (-1)^(xorder+1))/(xorder+1);
end

% matrix in front of d/dt
dt = zeros(N+1);

for j = 0:N

  Lj = x_coeffs(:,j+1);
  
  xj = x(j+1);
  x_p_vec=zeros(N+1,1);
  for i = 0:N
    x_p_vec(i+1) = xj^i;
  end
  Lj_xj = Lj.*x_p_vec; 
  % [x_p_vec Lj Lj_xj]
  
  integral = jac*wts(j+1)*sum(Lj_xj)*sum(Lj_xj);
  dt(j+1,j+1) = integral;

end

%dt

%% Analytically compute integral
dt2 = zeros(N+1);

for j = 0:N

Lj = x_coeffs(:,j+1);

integral = 0;
xj = x(j+1);
for i = 0:N
  for k = 0:N
    xorder=k+i;
    coeff = Lj(i+1)*Lj(k+1);
    integral = integral + coeff*dx_integral(xorder+1);
  end
end

dt2(j+1,j+1) = jac*integral;

end


% matrix for vdu/dx           % variational form: -udv/dx
dc = zeros(N+1);

for j = 0:N      % Test function (v)

  DLj = D(:,j+1);

  for i = 0:N      % Trial function (u)
  
    Li = x_coeffs(:,i+1);
    x_pi_vec=zeros(N+1,1);
    xi = x(i+1);
    for k = 0:N                   
        x_pi_vec(k+1) = xi^k;
    end
    x_pj_vec=zeros(N+1,1);
    for k = 1:N                   
        x_pj_vec(k+1) = xi^(k-1);
    end
  
    Li_xj = Li.*x_pi_vec;
    dLj_xj = DLj.*x_pj_vec; 
  
  %  if j == 2
  %       [x_p_vec dLi dLi_xj]
  %  end
    integral = jac*dxibydx*wts(i+1)*sum(Li_xj)*sum(dLj_xj);
  
    dc(j+1,i+1) = integral;
  
  end

end

dc = -dc;


%% Boundary term  for convection        % after integration by parts: uv(1) - uv(-1)
cbc = zeros(N+1);
loop = 1;


for j = 0:loop      % Test function (v)
     
     Lj = x_coeffs(:,j*N+1);

     xj = x(j*N+1);
     x_pj_vec=zeros(N+1,1);
     for k = 0:N
         x_pj_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_pj_vec; 
     % [x_p_vec Lj Lj_xj]

for i = 0:N      % Trial function (u)

     Li = x_coeffs(:,i+1);
     x_pi_vec=zeros(N+1,1);
     for k = 0:N 
         x_pi_vec(k+1) = xj^k;
     end
     Li_xj = Li.*x_pi_vec;

     integral = jac*((-1)^(j+1))*sum(Li_xj)*sum(Lj_xj);   % (-1)^(j+1) - cheap way to specify outgoing normal to surface.

     cbc(j*N+1,i+1) = integral;

end

end

cbc = cbc;
%

dc2 = zeros(N+1);                       % convection operator without integration by parts:    vdu/dx

for j = 0:N      % Test function (v)

     Lj = x_coeffs(:,j+1);

     xj = x(j+1);
     x_p_vec=zeros(N+1,1);
     for k = 0:N
         x_p_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_p_vec; 
     % [x_p_vec Lj Lj_xj]

     integral = 0;
for i = 0:N      % Trial function (u)

     dLi = D(:,i+1);
     x_p_vec=zeros(N+1,1);
     for k = 1:N                   % k=0 ==> has coefficient of 0. Hence skipped.
         x_p_vec(k+1) = xj^(k-1);
     end
     dLi_xj = dLi.*x_p_vec;
%     if j == 2
%          [x_p_vec dLi dLi_xj]
%     end
     integral = jac*dxibydx*wts(j+1)*sum(Lj_xj)*sum(dLi_xj);

     dc2(j+1,i+1) = integral;

end

end

dc2 = dc2;

%% Spatially inhomogenous convection
% matrix for vCdu/dx           % variational form: -udv/dx
VCdx = zeros(N+1);             % Just testing with uniform for now
xd_p_mat = zeros(Nxd+1,N+1);
xd_dp_mat = zeros(Nxd+1,N+1);


% Powers and power-1
for k=0:N
  xd_p_mat(:,k+1) = transpose(xd.^k);

  if (k==0)
    xd_dp_mat(:,k+1) = zeros(Nxd+1,1);
  else
    xd_dp_mat(:,k+1) = transpose(xd.^(k-1));
  end    
     
end
deriv_2 = transpose(transpose(D)*transpose(xd_dp_mat));
test_f  = transpose(transpose(x_coeffs)*transpose(xd_p_mat));
xgll = xs + (xe-xs)/2*(xd+1);
%Conv_xd = 1.0 + 0.1*sin(pi/2*xgll)
%Conv_xd = 1.0 + 0.1*xgll;
Conv_xd = ones(Nxd+1,1);

for j = 0:N    % Test function (v)

  for i = 0:N  % Trial function (u)
   
    sum_ij = Conv_xd.*test_f(:,j+1).*deriv_2(:,i+1).*wtsd;  

    integral = jac*sum(sum_ij);
  
    VCdx(j+1,i+1) = integral;
  
  end

end

VCdx = VCdx;


%% Laplacian term                  % Integration by parts:  dudx.dvdx
lp = zeros(N+1);

for j = 0:N                   % Test function (v)
     dLj = D(:,j+1);
for i = 0:N                   % Trial function (u)
     dLi = D(:,i+1);
     
     integral = 0;
for k = 0:N

     xk = x(k+1);

     x_p_vec=zeros(N+1,1);
     for kk = 1:N                   % k=0 ==> has coefficient of 0. Hence skipped.
         x_p_vec(kk+1) = xk^(kk-1);
     end

     dLi_xk = dLi.*x_p_vec;
     dLj_xk = dLj.*x_p_vec; 
     
     integral = integral + (jac*(dxibydx)^2)*wts(k+1)*sum(dLj_xk)*sum(dLi_xk);

end

     lp(j+1,i+1) = integral;

end

end

lp = -lp;


lpbc = zeros(N+1);             % Boundary terms for laplacian operator. vdudx(1) - vdudx(-1)
loop = 1;
for j = 0:loop      % Test function (v)
     
     Lj = x_coeffs(:,j*N+1);

     xj = x(j*N+1);
     x_pj_vec=zeros(N+1,1);
     for k = 0:N
         x_pj_vec(k+1) = xj^k;
     end
     Lj_xj = Lj.*x_pj_vec; 
     % [x_p_vec Lj Lj_xj]

for i = 0:N      % Trial function (u)

     dLi = D(:,i+1);
     x_pi_vec=zeros(N+1,1);
     for k = 1:N                   % derivative. hence k=0 skipped.
         x_pi_vec(k+1) = xj^(k-1);
     end
     dLi_xj = dLi.*x_pi_vec;

     integral = jac*dxibydx*((-1)^j)*sum(dLi_xj)*sum(Lj_xj);

     lpbc(j*N+1,i+1) = integral;

end

end

lpbc = lpbc;



% Hack for testing ...
%for i = 1:N+1
%     dc2(i,i) = 0;
%end

%%

%% Build forcing

%% Build basis change matrix
%% Taken directly from nek
nx = N+1;
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
spectonodal = transpose(pht);
nodaltospec = inv(spectonodal);
%utest = ones(N+1,1);
%leg_coeff = nodaltospec*utest;

phi = pht;

boyd = 0;
if boyd
  for i = 3:nx
    for j = 1:nx
      phi(i,j) = pht(i,j) - pht(i-2,j);
    end
  end
end
boydstonodal = transpose(phi);
nodaltoboyds = inv(boydstonodal);


% Filter function. Modes and amplitudes
G = zeros(N+1);
G(N+1,N+1)=0.25;
% Testing
%uleg_test = zeros(N+1,1);
%uleg_test(N+1,1)=1;
%unodal = spectonodal*uleg_test;
%
%boyds = G*nodaltoboyds*unodal;
%unodal2 = boydstonodal*boyds;
%uleg2 = nodaltospec*unodal2;



%%
% matrix for filter forcing

%G(N,N) = 0.5;
%G(N,N) =1;
%G(N-1,N-1) =0.1;

%G = eye(N+1);

% Build matrix which applies filter on Boyd transformed basis.
% And then transforms back to physical space.
fil_mat = boydstonodal*G*nodaltoboyds;

fil = zeros(N+1);

for j = 0:N      % Test function (v)

  Lj = x_coeffs(:,j+1);

  xj = x(j+1);
  x_p_vec=zeros(N+1,1);
  for k = 0:N
    x_p_vec(k+1) = xj^k;
  end
  Lj_xj = Lj.*x_p_vec; 

  Fi = fil_mat(j+1,:);
  
  integral = jac*wts(j+1)*sum(Lj_xj).*Fi;

  forc(j+1,:) = integral;

end


return
