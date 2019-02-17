% First attempts at coarsening the mesh

clear
clc
close all

%load saab_wing2d.mat
%load saab750k.mat
%load fluent_plus2.mat
%load stretched.mat
%load saab600k.mat
%load naca0012_5_15.mat
load naca0012_5_15_2.mat
%load naca0009_3p5c_c.mat

skiplayers = 2;         % Need to skip some layers since its smaller than the others
curvedef   = 'mv ';
ifvtk      = 1;
rea2d = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk);

ifre2 = 0;
Nek_WriteRea(rea2d,ifre2);



