function [ eigvec, M ] = setM( eigvec, W)
% [ eigvec, M ] = setM( eigvec, W)
% This routine sets up the matrix M_ij=<q_i,q_j>
%
% (c) Ardeshir Hanifi 2017


W4=blkdiag(W,W,W,zeros(size(W)));

% compute the inner product
E=diag(eigvec'*W4*eigvec);

% scale the eigenfunctions
eigvec=eigvec*diag(1./sqrt(E));

% recompute the inner product
M=eigvec'*W4*eigvec;
end

