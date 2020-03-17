function solnplot(El,nelv,h)

     figure(h)
     hold off
     for elno=1:nelv
          figure(h)
          surf(El(elno).xm1,El(elno).ym1,El(elno).soln);
          colorbar;
          hold on
     end
     hold off
