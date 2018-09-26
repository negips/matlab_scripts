%     First attempts at modifying the mesh

clear
clc
close all

casename = 'lu';

rea = Nek_ReadRea(casename);

n=rea.mesh.nelg;
cmap = jet(n);

figure(1);
for i=1:n
  fill(rea.mesh.xc(:,i),rea.mesh.yc(:,i),cmap(i,:)); hold on
end
colormap

