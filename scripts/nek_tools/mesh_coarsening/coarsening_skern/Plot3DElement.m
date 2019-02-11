function Plot3DElement(mesh3d,i,fig)

  figure(fig)
  
  ind=[1 2 3 4 1];
  plot3(mesh3d.xc(ind,i),mesh3d.yc(ind,i),mesh3d.zc(ind,i), 'k'); hold on

  ind=[5 6 7 8 5];
  plot3(mesh3d.xc(ind,i),mesh3d.yc(ind,i),mesh3d.zc(ind,i), 'k'); hold on

  ind=[1 5 6 2 1];
  plot3(mesh3d.xc(ind,i),mesh3d.yc(ind,i),mesh3d.zc(ind,i), 'k'); hold on

  ind=[4 8 7 3 4];
  plot3(mesh3d.xc(ind,i),mesh3d.yc(ind,i),mesh3d.zc(ind,i), 'k'); hold on

  xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16)
  ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16)
  zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 16)
  view(3)

end
%---------------------------------------------------------------------- 
