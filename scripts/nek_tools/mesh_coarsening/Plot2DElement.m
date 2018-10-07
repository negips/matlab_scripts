function Plot3DElement(mesh3d,i,fig)

  figure(fig)
  
  ind=[1 2 3 4 1];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  ind=[5 6 7 8 5];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  ind=[1 5 6 2 1];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  ind=[4 8 7 3 4];
  plot3(mesh3d.XC(ind,i),mesh3d.YC(ind,i),mesh3d.ZC(ind,i), 'k'); hold on

  xlabel('x')
  ylabel('y')
  zlabel('z')
  view(3)

end
%---------------------------------------------------------------------- 
