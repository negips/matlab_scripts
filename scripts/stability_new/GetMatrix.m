function EM = GetMatrix(OS,M,k2);
%
% This function computes the energy weighted 
% OS equation
%
% INPUT
% OS       = 3D Orr-Sommerfeld operator 
% M        = energy weight matrix
% k2       = alpha^2+beta^2
%
% OUTPUT 
% EM       = energy-weighted OS matrix 
 
    % Phase 1: Compute eigenvalues and eigenfunctions of 
    % Orr-Sommerfeld matrix and sort in order of descending 
    % imaginary part. The function nlize normalizes the 
    % eigenfunctions with respect to the weight matrix M.

    [xs,es] = OrderedEig(OS);
    xs      = ENormalize(xs,M);

    % Phase 2: Choose which modes are to be used in optimal 
    % calculation. Modes with imaginary part > 1 are neglected. 
    % Modes with imaginary part < imin are neglected as well.

%    dbstop in GetMatrix at 27   

    ishift = 1; 
    imin   = -1.5;

    while imag(es(ishift))>1, 
      ishift = ishift+1; 
    end

    [n1,n2] = ExtractEig(es,imin);

    cols    = (ishift:n2);
    xu      = xs(:,cols);
    eu      = es(cols); 
    ncols   = length(cols);
    fprintf('Number of modes used: %1.0f \n',ncols); 

    % Phase 3: Compute the reduced Orr-Sommerfeld operator

    [qb,invF] = TMatrix(M,xu,eu);
    EM        = sqrt(-1)*qb;
    
