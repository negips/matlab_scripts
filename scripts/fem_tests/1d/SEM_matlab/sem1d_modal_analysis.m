% Modal analysis of a 1D medium discretized by the Spectral Element Method
% Scalar wave equation
%
% INPUT:	X(:)	physical coordinates of the elements vertices (ordered)
%		RHO	object for the generation of density distribution
%			RHO.fun(RHO.data,e,x) gives the density 
%			on the GLL points x inside element e
%			See examples in sem1d_homog.m and sem1d_two_layers.m
%		MU	same for the shear modulus
%		NGLL	number of Gauss-Lobatto-Legendre points per element 
%			= polynomial degree +1
%		BC	['NN'] boundary conditions at min(X) and max(X) respectively
%				N	homogeneous Neumann (stress free)
%				D	homogeneous Dirichlet (no displacement)
%
function [eigvec,eigval,x,err]=sem1d_modal_analysis(X,RHO,MU,NGLL,BC)

%---- Problem setup ----

% Read the Gauss-Lobatto-Legendre nodes, quadrature weights
% and derivatives of lagrange interpolants Hij=hprime_i(xi_j)
[xgll,wgll,H] = GetGLL(NGLL);

nel = length(X)-1;	% total number of elements
nglob = nel*(NGLL-1)+1; % total number of nodes

% initialize vectors and matrices
K = zeros(nglob,nglob,1);	% stiffness matrix, 
				% Note: for large scale problems it is smarter 
				% to store K as sparse
M = zeros(nglob,1);		% mass matrix
x = zeros(nglob,1);		% coordinates of the computational nodes

% build x, M and K :
for e=1:nel, 		% loop over elements

 % global indices associated to the GLL nodes of this element:
  iglob=(e-1)*(NGLL-1)+(1:NGLL); 

 % local to global coordinate map (xi --> x)
  x(iglob) = 0.5*(1-xgll)*X(e) + 0.5*(1+xgll)*X(e+1); 
  dxe=0.5*(X(e+1)-X(e));	% dx/dxi

 % for the sake of generality, the physical properties are set
 % by user-defined functions and data
  mu  = feval(MU.fun,MU.data,e,x(iglob));
  rho = feval(RHO.fun,RHO.data,e,x(iglob)); 	

 % mass matrix assembly
 % local contribution = rho(:)*GLL_weights(:)*dx/dxi
  M(iglob) = M(iglob) + rho .*wgll *dxe;

 % stiffness matrix assembly
 % local contribution = H*(mu(:)*GLL_weights(:)*dxi/dx)*H_transposed
  K(iglob,iglob) = K(iglob,iglob) + H * ( repmat(mu.*wgll/dxe,1,NGLL).* H');

end



%---- Boundary Conditions ----

% Default: homogeneous Neumann (stress free) on both ends
if ~exist('BC','var'), BC='NN'; end 

% modify matrices if any Dirichlet condition: 
% eliminate the boundary degree of freedom
if ~strcmp('NN',BC)
  if BC(1)=='D', i1=2; else, i1=1; end
  if BC(2)=='D', i2=nglob-1; else, i1=nglob;  end
  M = M(i1:i2);
  K = K(i1:i2,i1:i2);
  nglob = i2-i1+1;
end




%---- Eigenvalue solver ----

% We're set up now to compute the eigenvalues and eigenvectors of M^(-1)*K 
% However, the eigenvalue solver works best with symmetric matrices
% and M^(-1)*K is NOT symmetric.
% The trick is to define a new symmetric matrix:   A = M^(-1/2)*K*M^(-1/2)
% The eigenvalues of the original problem are the same as those of A,
% the eigenvectors are V = M^(-1/2) * V_of_A

% building A = M^(-1/2)*K*M^(-1/2)
invsqrtM = 1./sqrt(M);
A = K .* repmat(invsqrtM',nglob,1);
A = repmat(invsqrtM,1,nglob) .*A;
A = 0.5*(A+A'); % just to be sure A is exactly symmetric

% Find all eigenvalues/eigenvectors.
% The standard Matlab routine 'eig' does not exploit sparseness :(
% The routine 'eigs' does, but only if a few eigenpairs are requested
[eigvec,eigval]=eig(A);

eigval=sqrt(abs(diag(eigval))); 	% eigenvalues of A are frequency^2

eigvec=repmat(invsqrtM,1,nglob) .* eigvec;	% rescale eigenvectors

% output in ascending eigenvalue order
[eigval,isort]=sort(eigval);
eigvec= eigvec(:,isort);

% if any Dirichlet condition: add the boundary node (=0) to eigenvectors
if ~strcmp('NN',BC)
  pad = zeros(1,nglob);
  if BC(1)=='D', eigvec = [ pad; eigvec ]; end 
  if BC(2)=='D', eigvec = [ eigvec; pad ]; end
end

% Very low frequency eigenfrequencies are not well recovered
% due to round-off error.
% This is an estimate of the absolute error on eigenvalues = err./eigval
err=0.5*eps*norm(A,2);
