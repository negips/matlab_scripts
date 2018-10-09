function Plot3DElement(mesh2d,i,fig)

  figure(fig)
  
  ind=[1 2 3 4 1];
  plot(mesh2d.xc(ind,i),mesh2d.yc(ind,i), 'k'); hold on

  xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16)
  ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16)

end
%---------------------------------------------------------------------- 
