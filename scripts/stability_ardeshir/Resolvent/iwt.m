function W = iwt(N)
% IWT: sets up the weight function for the spectral integration and
% form the matrix used for computation of disturbance energy.
%
% Input
%    N: order of Chebychev olynomial is N-1
%
% Output:
%    W: array containing the integration weights
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%


N1=N-1;
%compute chebycheff integration weights
n = 0:1:N1;
j = 0:1:N1;
b = ones(1, N);
b([1 N]) = 0.5;
c = 2*b;
b = b/N1;
S = cos(n(3:N)'*j*(pi/N1));
W = diag(b.*((2+(c(3:N).*((1+(-1).^n(3:N))./(1-n(3:N).^2)))*S)));

end
