% 3D mesh generation

clear
clc
close all

load saab_wing2d.mat
%load saab750k.mat
%load fluent_plus2.mat
%load stretched.mat

skiplayers = 3;         % Need to skip some layers since its smaller than the others
curvedef   = 'mv ';     % Definition of curved boundaries
ifvtk      = 1;
rea2d = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk);


% Generate 3D mesh
nz0=8;
Lz=0.15;
ifperiodic=1;
ifvtk=1;
rea3d = Generate3DCoarse(rea2d,nz0,Lz,ifperiodic,ifvtk); 

Nek_WriteRea(rea3d,1);

