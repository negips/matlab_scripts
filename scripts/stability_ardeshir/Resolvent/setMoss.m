function [ eigvec, M ] = setMoss( eigvec, W, D, alpha, beta)
% [ eigvec, M ] = setMoss( eigvec, W, D, alpha, beta)
% This routine sets up the matrix M_ij=<q_i,q_j>
%
% (c) Ardeshir Hanifi 2017

D2=D(:,:,2);

I=eye(size(W));
W2=blkdiag(W,W);
k2=alpha^2+beta^2;

% compute the inner product
F=W2*blkdiag(I-D2/k2,I/k2);
E=diag(eigvec'*F*eigvec);

% scale the eigenfunctions
eigvec=eigvec*diag(1./sqrt(E));

% recompute the inner product
M=eigvec'*F*eigvec;

end

