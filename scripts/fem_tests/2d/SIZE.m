% Declarations for the case
% Trying to follow Nek.
% Lots of things are different of course.

% Polynomial orders

Nx = 12;
lx1=Nx+1;
%Nxd = 9;
Nxd = ceil(3/2*(Nx+1));

Ny = Nx;            % For now. Maybe will test a more general case later.
ly1=Ny+1;
Nyd = Nxd;

display(['Nx=' num2str(Nx) '; Nxd=' num2str(Nxd)])
display(['Ny=' num2str(Ny) '; Nyd=' num2str(Nyd)])

