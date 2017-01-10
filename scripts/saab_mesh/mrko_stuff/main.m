%Generate geometry for the nlf wing
clear all
close all
addpath ./wing_tools/
addpath('/scratch/mrko/work/matlab-library/nek')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Free-stream velocity: 30 m/s
%      Angle of attack     : -5 degrees
%      Sweep angle         : 45 degrees
%      Chord length        : 0.7 m (normal to leading edge).

loadgrid = true;
%0)UVP Import file    (check in gen_import_uvp)
importfile8 = 'wing25D0.f00027'; %double precision
    % importfile4 = 'wing25D0.f00112'; %single precision
    %INTERPOLATION also possible in nek with g2gi
    % importfile_g2gi = 'new0.f00001'; %single precision in /g2gi2D/
%elements in spanwise direction
nelz = 5;  %depending on extrusion through n2to3

% 1) Geometry
Dim = 2;
wingdata = './buterfli/wing.mat';
sweep    = 60;
short_chord_l = 0.35;

% qref = 70;
% -- exponent for weight function
rexp = 0.2;%.2;            %relating dx to local prof curv

% -- number of points for the grid that will be mapped
nelx = 200;%200 100           % mapping grid along the profile
nely =  18;%50 18             % parameter for elements in y-dir./rectgen  


% ReC = 4.399864912605778e+06;
% aoa = 0;                  % Angle of attack
% swa = 0.;                  % Sweep angle
% Cchord =1;%1.83*cos(swa);

frac_chord_up = 50;%50        % Portion of wing ([%] of short chord), upper
frac_chord_lo = 1.0;%1.5;     % Portion of wing ([%] of chord), lower
hout = .015;%.02                   %outflow height
mult          = 5;            % Refinement of airfoil cont.; 
                              % <mult> splittings between given pts.
% chord         = 1.0;

refmethod = 2;              % Refinement method (tang. dir.): 1 - cos; 2 - geometric
aluratio = 1.1;

refsystem = 1;              % Ref. length: 1 - short chord; 2 - nose radius

% % 2) SEM-mesh parameters
% nelxuref = 40;              % # elements in refined nose region
% frac_chord_up_fine  = .15;  % Refined portion of wing (relative to frac_chord_up*chord)
% frac_chord_low_fine = .15;

% 3) setting for the element numbers and orders
nelx_scale_elnum = 0.45;
nely_scale_elnum = 0.3;
% nely = 35;  
% N = 6;                     % Polynomial order
             

% mesh qulity variabl
mesh_satisfactory = false;
run_mesh = 1;
% nekton input files
casename ='wing2Dsmall';
casename3D = 'wing3Dsmall';
casename3D_g2gi = 'wing3Dsmallg2gi';

path = './';
pathsmall = './small/';

color = ['b','g','m','r','y','k','b','m','g','r','k','y'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get outer boundary based on streamlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define box dimensions
curve_in = false;

xl = -0.5;          %x-coordinate of left straight boundary
yl1 = -0.1;%-0.2;   %-5.6904;      %y-coordinate of lower left corner
yl2 = 0.2;%20;       %y-coordinate of upper left corner

xl3 = -0.1;
yl3 = 0;

xr1 = -5.2;  %9.2;          %x-coordinate of lower outlet boundary
xr2 = 10.1;         %x-coordinate of upper outlet boundary
  

ystr1 = -.12;%-5.7315;   %-5.6904;      %y-coordinate of lower left corner
ystr2 = 0.2;       %y-coordinate of upper left corner


%%
if ~loadgrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% get geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xy_coarse,xy_fine,m_coarse,nu_coarse,nl_coarse,nelx,lmult,m0] = ...
%     getgeom(aoa,swa,Cchord,frac_chord_up,frac_chord_lo,mult,chord,...
%     refmethod,aluratio,refsystem,frac_chord_up_fine ...
%     ,frac_chord_low_fine,nelxuref); 

%% Load original profile coordinates

data = load(wingdata); 
xpro = data.wing(:,1)'*cosd(sweep); 
shift = min(xpro);
xpro = xpro-shift;
ypro = data.wing(:,2)';
% xpro = fliplr(xpro);
% ypro = fliplr(ypro);
figure;plot(data.wing(:,1)-min(data.wing(:,1)),data.wing(:,2),'k',xpro,ypro,'r*-')
axis equal; legend('long chord profile','short chord profile')

% tangential coordinate
dx = diff(xpro); 
dy = diff(ypro); 

ds = sqrt(dx.^2+dy.^2);
spro = cumsum([0 ds]); 
[~,ispro0] = min(xpro); 
spro = spro - spro(ispro0);

% profile spline
xprfun = csapi(spro,xpro);
yprfun = csapi(spro,ypro);

% downsample
spro = linspace(spro(1),spro(end),100000);
xpro = fnval(xprfun,spro);
ypro = fnval(yprfun,spro);

xprfun = csapi(spro,xpro);
yprfun = csapi(spro,ypro);

% curvature radius (spline differentiation)
d1x = fnval(fnder(xprfun,1),spro); 
d2x = fnval(fnder(xprfun,2),spro);
d1y = fnval(fnder(yprfun,1),spro); 
d2y = fnval(fnder(yprfun,2),spro);

rpro = abs( (d1x.^2 + d1y.^2).^1.5 ./ (d1x.*d2y - d2x.*d1y)   );

rprfun = csapi(spro,rpro);

% tangential and normal directions
tpro = [  d1x ; d1y ]; 
tpro = tpro ./ ([1;1]*(sqrt(tpro(1,:).^2 + tpro(2,:).^2)));
npro = [ -d1y ; d1x ]; 
npro = npro ./ ([1;1]*(sqrt(npro(1,:).^2 + npro(2,:).^2)));
      
tprfun{1} = csapi(spro,tpro(1,:)); 
tprfun{2} = csapi(spro,tpro(2,:));
nprfun{1} = csapi(spro,npro(1,:)); 
nprfun{2} = csapi(spro,npro(2,:));

%% Cut and re-interpolate profile
xprcutlw = frac_chord_lo*short_chord_l*0.01;
xprcutup = frac_chord_up*short_chord_l*0.01;
% cut profile
% icb  = find(xpro < xprcutlw,1,'first');
icb  = find(xpro < xprcutlw,1,'last');
ice  = find(xpro < xprcutup,1,'last'); 
scut = spro([icb,ice]);

scut(1) = fzero(@(s) fnval(xprfun,s)-xprcutlw,scut(1));             %check
scut(2) = fzero(@(s) fnval(xprfun,s)-xprcutup,scut(2));

% new spacing: weight function (wpr) for ds based on the profile curvature (1/rpr)
wprfun = fnint(csapi(spro,abs(1./smooth(smooth(rpro))).^(rexp)));
wpro   = fnval(wprfun,spro);
wprinv = csapi(wpro,spro);

wpr = linspace(fnval(wprfun,scut(1)),fnval(wprfun,scut(2)),nelx+1);

% interpolation
spr = fnval(wprinv,wpr);
xpr = fnval(xprfun,spr);
ypr = fnval(yprfun,spr);
tpr =[fnval(tprfun{1},spr);
      fnval(tprfun{2},spr)];
npr =[fnval(nprfun{1},spr);
      fnval(nprfun{2},spr)];
rpr = fnval(rprfun,spr);

ipr = 1 + 0*xpr;

figure;plot(xpro,ypro,'k.',xpr,ypr,'r*');axis equal

xy_coarse(:,1) = xpr;
xy_coarse(:,2) = ypr;
xy_fine(:,1) = xpro;
xy_fine(:,2) = ypro;

%
%%
%grid measures
% - upper left corner of grid / start
 xs = xpr(1) + hout*npr(1,1);
 ys = ypr(1) + hout*npr(2,1);
 
% - upper right corner of grid /end
 xe = xpr(end) + hout*npr(1,end);
 ye = ypr(end) + hout*npr(2,end);

 xo    = linspace(xs,xe,max(size(xpr)));%nxul);
    
 xous = (xpr-min(xpr));
 xous = xous./max(xous);
 xous = xous.*(xe+abs(xs))+xs;
 
 yous = (ypr-min(ypr));
 yous = yous./max(yous);
 yous = yous.*(ye-ys)+ys;

 yo = interp1(xous,yous,xo);

%end
hold on
plot(xo,yo,'b*-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute domain height at lower and upper outlet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
height_in = sqrt((xo(1)-xy_coarse(1,1))^2 + (yo(1) - xy_coarse(1,2))^2);
height_out = sqrt((xo(end)-xy_coarse(end,1))^2 + (yo(end) - xy_coarse(end,2))^2);

disp(['Domain height of lower outlet: ' num2str(height_in)])
disp(['Domain height of upper outlet: ' num2str(height_out)])

%--------------------------------------------------------------------------
%% plot whole geometry
%--------------------------------------------------------------------------
npx = nelx+1;
npy = nely+1;

xg = [xy_coarse(:,1)' fliplr(xo)];
yg = [xy_coarse(:,2)' fliplr(yo)];

%GE defines the outer boundary of the mesh
GE = [xg' yg'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create rectangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [GR,npxg,npyg] = rectangle_gen_small(nelx_scale_elnum,nely_scale_elnum,npx,npy);
% % [GR,npxg,npyg] = rectangle_gen_small62x23(nelx_scale_elnum,nely_scale_elnum,npx,npy);
% [GR,npxg,npyg] = rectangle_gen_small204x23y(nelx_scale_elnum,nely_scale_elnum,npx,npy);
[GR,npxg,npyg] = rectangle_gen_smallfocle(nelx_scale_elnum,nely_scale_elnum,npx,npy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% turnpoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_turnpoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Executing gridgen with the current settings + Visualization of the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system('gridgen prm.nlf')

%viewgrid('prm.nlf','v')
[x,y] = viewgrid('prm.nlf','v');
%vg('grid.nlf')
xp = reshape(x,npyg,npxg);
yp = reshape(y,npyg,npxg);
%%
% for i=1:size(xp,1)
%     plot(xp(i,:),yp(i,:),'r')
%     hold on
%     pause(.1)
% end

%% commented out for cut corners ***TEST***
% ind_upcorner = find(diff(xp(end,:))==0,1,'first');
% % ind_upcorner = find(diff(yp(end,:))<0,1,'first');
% ind_locorner = find(diff(xp(end,:))==0,1,'last' )+1;
% yp(end,1:ind_upcorner) = yup;
% xp(end,ind_upcorner) = xs;
% yp(end,ind_locorner:end) = ylow;
% xp(end,ind_locorner) = xs;

%%
% %Figure 99 is used to visualize the RANS simulation, the grid in this
% %plot is also plotted on top of the RANS simulation.
% fig6=figure(99);
% hold on
% for i = 1:npyg
%     plot(xp(i,:),yp(i,:),'r')
% end
% for i = 1:npxg
%     plot(xp(:,i),yp(:,i),'r')
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check and modify the mesh quality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh_check


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correct positions of vertices at airfoil boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Correct wall points')
checkobj = zeros(1,npxg);
for i = 2:npxg-1
    [box(1).x(1,i),box(1).y(1,i),checkobj(i)] = findbound2(box(1).x,box(1).y,xy_fine,i);
end
figure(200)
plot(checkobj)



%--------------------------------------------------------------------------
%% correct positions of the upper and lower bounds for ON boundary condition
%--------------------------------------------------------------------------
% fix4on_v2

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extend box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if extend_box
%    disp('Adding additional boxes')
%    box_extender_v2
% end

%--------------------------------------------------------------------------
%flip arrays such that first point is located at the lower part of the wing
%--------------------------------------------------------------------------
for n = 1:length(box)
  box(n).x = fliplr(box(n).x);
  box(n).y = fliplr(box(n).y);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
create_ELsmall
plot_geom


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the boundary condtions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
define_bcsmall


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixtail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixtail
% plot_geom

save grid
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write rea-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Dim == 2
%     writerea(EL,ElmConn,Nelm,path,FluidBC)
%     writerea2D_nocurve(EL,ElmConn,Nelm,path,FluidBC)
    writerea2D_nocurvesmall(ELsmall,ElmConnsmall,Nelm,path,FluidBC)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    load grid
    plot_geom
%     if Dim == 2
% %         writerea2D(EL,ElmConn,Nelm,path,FluidBC)
% writerea2D_nocurve(EL,ElmConn,Nelm,path,FluidBC)
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate and Import the gll points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gen_import_gllsmall


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extrude from small and base file (big domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeic_2D_small
% writeic_2D_small_inh
% extrude_smallto3D
extrude_smallto3D_T
% extrude_smallto3D_T_inh

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extrude from 2D g2gi file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extrude_smallto3D_g2gi


end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define fringe regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write_lambda2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%creat folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case_folder = ['wing' '_EL_' num2str(ceil(length(EL)/1000)) 'k'];
% display(['moving items to folder: ' case_folder ' ...']);
% system(['mkdir ' case_folder]);
% system(['cp ' casename '.* ' case_folder]);
% system(['cp mesh.rea ' case_folder]);
% system(['cp grid.mat ' case_folder]);
% display('Done!')
