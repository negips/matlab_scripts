function PlotVar(El,varname,nelv)

  h=figure;
  for elno=1:nelv
    evalstr = ['El(' num2str(elno) ').scrtch1 = El(' num2str(elno) ').' varname ';'];
    eval(evalstr)

    figure(h)
    surf(El(elno).xm1,El(elno).ym1,El(elno).scrtch1, 'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on
  end
  colorbar



