% 3D mesh generation

clear
clc
close all

%load saab_wing2d.mat
%load saab750k.mat
%load fluent_plus2.mat
%load stretched.mat
%load saab600k.mat
%load naca0012_5_15.mat
%load naca0012_5_15_2.mat
load naca0009_3p5c_c.mat

skiplayers = 2;         % Need to skip some layers since its smaller than the others
curvedef   = 'w  ';     % Definition of curved boundaries
ifvtk      = 1;
rea2d = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk);

%ifre2=0;
%Nek_WriteRea(rea2d,ifre2);

%return
% Generate 3D mesh
nz0=20;
Lz=0.10;
ifperiodic=1;
ifvtk=1;
rea3d = Generate3DCoarse(rea2d,nz0,Lz,ifperiodic,ifvtk); 

%ifre2 = 1;
%Nek_WriteRea(rea3d,ifre2);


% testing
%new=Nek_ReadRea('new');
