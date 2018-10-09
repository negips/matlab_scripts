% 3D mesh generation

clear
clc
close all

%load saab_wing2d.mat
%load saab750k.mat
load fluent_plus2.mat

skiplayers = 2;         % Need to skip some layers since its smaller than the others
curvedef   = 'W  ';     % Definition of curved boundaries
ifvtk      = 1;
mesh2d = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk);


% Generate 3D mesh
nz0=4;
Lz=1.0;
ifperiodic=1;
ifvtk=1;
mesh3d = Generate3DCoarse(mesh2d,nz0,Lz,ifperiodic,ifvtk); 


