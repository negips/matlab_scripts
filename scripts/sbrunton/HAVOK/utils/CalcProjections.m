function [Xi resid]= CalcProjections(Theta,dXdt)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% compute Sparse regression: sequential least squares

% Use Gram-Schmid to project out modes.

[m n] = size(Theta);
resid = dXdt;
for ind = 1:n                   % n is state dimension
  Xi(ind) = sum(resid.*Theta(:,ind));
  resid = resid - Xi(ind)*Theta(:,ind);        
end
