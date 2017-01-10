% Testing 2D mass integration

clear
clc

Nx=6;
Ny=Nx;

%xc = [-1 1 1.2 -1.1];
%yc = [-1 -0.9 1.1 1];

xc = [-1 1 1 -1];
yc = [-1 -1 1 1];


mass = lastmodeint(Nx,Ny,xc,yc);
[S2N N2S] = Leg2Nodal(Nx);

Nxd=Nx+1;
Nyd=Nxd; 
massd = lastmodeintd(Nx,Ny,Nxd,Nyd,xc,yc);

spectra_x = zeros(Nx+1,1);
spectra_x(end) = 1;
spectra_y = zeros(1,Ny+1);
spectra_y(end) = 1;
nodal_x = S2N*spectra_x;
nodal_y = spectra_y*transpose(S2N);

spectra_tens = kron(spectra_x,spectra_y);
nodal_tens = kron(nodal_x,nodal_y);

inner_mass = mass*nodal_tens(:);
inner_massd = massd*nodal_tens(:);

mass_fld = reshape(inner_mass,Nx+1,Ny+1);
massd_fld = reshape(inner_massd,Nx+1,Ny+1);

figure
surf(mass_fld); colorbar;
title('Mass matrix')

figure
surf(massd_fld); colorbar;
title('Dealiased mass matrix')

ratios = mass_fld./massd_fld;
figure
surf(ratios); colorbar;
title('Nodal value ratios')


