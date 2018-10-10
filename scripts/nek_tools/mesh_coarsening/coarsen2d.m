% First attempts at coarsening the mesh

clear
clc
close all

%load saab_wing2d.mat
%load saab750k.mat
%load fluent_plus2.mat
load stretched.mat

skiplayers = 3;         % Need to skip some layers since its smaller than the others
curvedef   = 'mv ';
ifvtk      = 1;
rea2d = Generate2DCoarse(rea,LayerE,LayerX,LayerY,LayerBC,LayerCEl,MeshC,skiplayers,curvedef,ifvtk);


Nek_WriteRea(rea2d);

