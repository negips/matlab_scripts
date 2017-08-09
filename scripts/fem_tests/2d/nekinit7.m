% Small test simulation.

clear
clc
close all

%addpath '../../'

tstart = tic;

display(['Time started: ' datestr(clock)])

SIZE           %    Polynomial order
re2            %    domain/decomposition/mapping/boundary conditions.

Re=1e+5;
nu=1/Re;
ifboyd = 0;
destn  = 'plots/';
ifplot = 0;
ifdiv  = 1;
gltot  = 0;     %    Global number of elements 

for ii=1:nelv

     % clc 
     display(['Building Matrices for Element ', num2str(ii)])

     [mass nek_mass DXM1 DYM1 DXM1D DYM1D RXM1 RYM1 SXM1 SYM1 convx convy convall convxd convyd convalld convxd_new convyd_new convalld_new Cfx Cfy gradm1x gradm1y gradm1xd gradm1yd intpm1d wtsvecd nek_conv lpx lpy lpall nek_lp lpbc forc NxNy_nodal2spec NxdNyd_spec2nodal x_coeff y_coeff Dx Dy w2m1 xm1 ym1 JACM1 JACM1D xm1d ym1d] = MESem2D7(Nx,Ny,Nxd,Nyd,El(ii).xc,El(ii).yc,ifboyd,ifplot);

     El(ii).mass = mass;      
     El(ii).nek_mass = nek_mass;
     El(ii).Dx = Dx;
     El(ii).Dy = Dy;
     El(ii).x_coeff = x_coeff;
     El(ii).y_coeff = y_coeff;
     El(ii).DXM1 = DXM1;
     El(ii).DYM1 = DYM1;
     El(ii).DXM1D = DXM1D;
     El(ii).DYM1D = DYM1D;
     El(ii).RXM1 = RXM1;
     El(ii).RYM1 = RYM1;
     El(ii).SXM1 = SXM1;
     El(ii).SYM1 = SYM1;
     El(ii).convx = convx;
     El(ii).convy = convy;
     El(ii).convall = convall;
     El(ii).forc = forc;
     El(ii).NxNy_nodal2spec = NxNy_nodal2spec;
     El(ii).NxdNyd_spec2nodal = NxdNyd_spec2nodal;

     El(ii).convxd = convxd;
     El(ii).convyd = convyd;
     El(ii).convalld = convalld;

     El(ii).convxd_new = convxd_new;
     El(ii).convyd_new = convyd_new;
     El(ii).convalld_new = convalld_new;

     El(ii).Cfx = Cfx;
     El(ii).Cfy = Cfy;

     El(ii).gradm1x = gradm1x;
     El(ii).gradm1y = gradm1y;
     El(ii).gradm1xd = gradm1xd;
     El(ii).gradm1yd = gradm1yd;
     El(ii).intpm1d = intpm1d;
     El(ii).wtsvecd = wtsvecd;
     El(ii).xm1d = xm1d;
     El(ii).ym1d = ym1d;
     El(ii).JACM1D = JACM1D;

     El(ii).nek_conv = nek_conv;
     El(ii).laplx = lpx;
     El(ii).laply = lpy;
     El(ii).lapl = lpall;
     El(ii).nek_lp = nek_lp;
     El(ii).xm1 = xm1;
     El(ii).ym1 = ym1;
     El(ii).JACM1 = JACM1;
     El(ii).w2m1 = w2m1;

     El(ii).un = 0*xm1;
     El(ii).ulag1 = El(ii).un;
     El(ii).ulag2 = El(ii).un;
     El(ii).soln = El(ii).un;
     El(ii).b = El(ii).un;
     El(ii).bl = El(ii).un;
     El(ii).scrtch1 = El(ii).un;
     El(ii).scrtch2 = El(ii).un;

     %% Calculate divergence
     div = El(ii).gradm1x*Cfx(:) + El(ii).gradm1y*Cfy(:);
     El(ii).div = reshape(div,Nx+1,Ny+1);
     divd = El(ii).gradm1xd*Cfx(:) + El(ii).gradm1yd*Cfy(:);
     El(ii).divd = reshape(divd,Nxd+1,Nyd+1);

     %%
end


%    Build multiplicity
for ii = 1:nelv
     [El(ii).lglmap gltot] = GlobalNo(ii,Nx,Ny,nelx,nely,xperiodic,yperiodic,gltot,El);
end
El = vmult(El,nelv,lx1,ly1);

%    inidices/Global nos of boundary nodes
El = BdryIndex(El,nelv,lx1,ly1);

% initial Condition
El = SetICs(El,nelv);

% Check Domain volume
volume = GetVolume(El,nelv);
disp(['Domain Volume: ' num2str(volume)])

% Check Divergence
divergence = GetDivNorm(El,nelv);
dnorm = divergence/volume;
disp(['Divergence Norm: ' num2str(dnorm)])

if ifdiv

  for ii=1:nelv

%   Check Divergence
    figure(100) 
    surf(El(ii).xm1,El(ii).ym1,log10(abs(El(ii).div)+1e-17), 'EdgeColor', 'none', 'FaceColor', 'interp'); hold on
    title('Divergence field (logscale)')
    view(2)
    colorbar  

%%   2D spectra of the convecting field. 
%    figure(101)
%    spectra2d = log10(abs(El(ii).NxNy_nodal2spec*El(ii).Cfx(:))+1e-10);
%    spectra2d = reshape(spectra2d,Nx+1,Ny+1);
%    surf(El(ii).xm1,El(ii).ym1,spectra2d, 'EdgeColor', 'none', 'FaceColor', 'interp'); hold on
%    title('Convecting field Cfx Spectra')
%    colorbar
%    view(2)
%
%%   2D spectra for the filtered Convecting field.
%    Just a check to see forcing matrix is correct.  
%    figure(102)
%    filtered_fld =  inv(El(ii).mass)*El(ii).forc*El(ii).Cfx(:);
%    spectra2d = log10(abs(El(ii).NxNy_nodal2spec*filtered_fld)+1e-20);
%    spectra2d = reshape(spectra2d,Nx+1,Ny+1);
%    surf(El(ii).xm1,El(ii).ym1,spectra2d, 'EdgeColor', 'none', 'FaceColor', 'interp'); hold on
%    title('Spectra filtered convective field')
%    colorbar
%    view(2)

%   Visualize Convecting field in x. 
%    figure(103)  
%%    dCdy = El(ii).gradm1yd*El(ii).Cfy(:);
%%    dCdy = reshape(dCdy, size(El(ii).xm1d));
%    convect = El(ii).Cfy;  
%    surf(El(ii).xm1,El(ii).ym1,convect, 'EdgeColor', 'none', 'FaceColor', 'interp'); hold on
%    title('Cy')
%    colorbar
%
%%   Visualize Convecting field in y 
%    figure(104)  
%    convection = El(ii).Cfy;
%    surf(El(ii).xm1,El(ii).ym1,convection, 'EdgeColor', 'none', 'FaceColor', 'interp'); hold on
%    title('Cfy')
%    colorbar

  end 

end  


plotgll=0;
%nekchecks(El,Nx,Ny,nelx,nely,nelv,plotgll)

plotspy=1;
[bigmass bigconv bigconvd bigconvxd bigforc bigconvd_new biglapl velvec gno nreps El] = AssembleBig7(El,Nx,Ny,nelx,nely,nelv,plotspy);

%% Check eigenvalues of the system
eigfigure=[];
ifplot=1;
if (ifplot)
  eigfigure=figure;
  hold on
end


sparsehandle=[];
ifsparse=1;
if (ifsparse)
  sparsehandle=figure;
  hold on
end

bdfkstability = 0;

%% Convective matrix with overintegration. eigen values  
sysmat = inv(bigmass +nu*biglapl)*bigconvd_new;
col='b';
[evec lambda] = SystemEig(sysmat,ifplot,eigfigure,ifsparse,sparsehandle,bdfkstability,col);
pause(2)    
%
%if (ifplot)
%  filename=['spectra_conv_N' num2str(Nx), '_Nxd' num2str(Nxd) '_nelv' num2str(nelv) '.eps'];
%  SaveFig(eigfigure,filename,destn,1);
%end
%
%%% (Convective - Forcing) matrix eigenvalues
%sysmat = inv(bigmass)*(bigconvd_new - bigforc);
%col='r';
%ifsparse=0;
%ifplot=1;
%[evec lambda] = SystemEig(sysmat,ifplot,eigfigure,ifsparse,sparsehandle,bdfkstability, col);
%pause(2)    
%
%if (ifplot)
%  filename=['spectra_rhs_N' num2str(Nx), '_Nxd' num2str(Nxd) '_nelv' num2str(nelv) '.eps'];
%%  SaveFig(eigfigure,filename,destn,1);
%end

%% Testing dealiased convection operator.
% Convective matrix eigen values 
sysmat = inv(bigmass)*bigconvd;
col='k';
ifsparse=0;
[evec lambda] = SystemEig(sysmat,ifplot,eigfigure,ifsparse,sparsehandle,bdfkstability,col);
pause(2)    
%
%%% (Convective - Forcing) matrix eigenvalues
sysmat = inv(bigmass + nu*biglapl)*(bigconvd - 0*bigforc);
col='r';
ifsparse=0;
ifplot=1;
[evec lambda] = SystemEig(sysmat,ifplot,eigfigure,ifsparse,sparsehandle,bdfkstability, col);
pause(2)    

%%

clearvars mass nek_mass DXM1 DYM1 DXM1D DYM1D RXM1 RYM1 SXM1 SYM1 convx convy convall convxd convyd convalld gradm1xd gradm1yd intpm1d wtsvecd nek_conv lpx lpy lpall nek_lp lpbc forc x_coeff y_coeff Dx Dy w2m1 xm1 ym1 JACM1 JACM1D xm1d ym1d

clearvars xpt ypt ee ii jj ifplot plotspy plotgll

display(['Initialization completed: ' datestr(clock)])
toc(tstart)

%save(['matrices_N' num2str(Nx) '_NELV' num2str(nelv) '_CART5_1.mat'], '-v7.3');
%% SOLVE
