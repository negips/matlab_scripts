% Renumber elements


coarsen3d

n2d=rea2d.mesh.nelg;
n3d=rea3d.mesh.nelg;
nz=n3d/n2d;

e=0;
for i=1:n2d
   for j=1:nz
      e=e+1;
      oldgl(e)=rea3d.mesh.globalno(e);
      newgl(e)=(j-1)*n2d+i;
   end
end

new=rea3d;

new.mesh.xc(:,newgl)=rea3d.mesh.xc;
new.mesh.yc(:,newgl)=rea3d.mesh.yc;
new.mesh.zc(:,newgl)=rea3d.mesh.zc;

new.mesh.cbc(:,newgl) = rea3d.mesh.cbc;

nfaces=6;
for i=1:n3d
   for j=1:nfaces
      c2=new.mesh.cbc(j,i).connectsto;
      if c2>0
        cnew=newgl(c2);
        new.mesh.cbc(j,i).connectsto=cnew;
      end  
   end
end

ncurves=new.mesh.Ncurve;
cieg = [];
for i=1:ncurves
   c2=new.mesh.curveieg(i);
   cnew=newgl(c2);
   cieg = [cieg cnew];
end   

[cieg2 ind]=sort(cieg);
cedge = new.mesh.curveedge(ind);

new.mesh.curveieg = cieg2;
new.mesh.curveedge= cedge;

CheckConnectivity3D(new.mesh)
%Nek_WriteRea(new,0)
     
