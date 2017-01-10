function Elout = UpdSoln(El,soln,gno,nelv,lx1,ly1,h,ifplot)

     figure(h)
     for elno=1:nelv
          El(elno).unvec = 0*El(elno).unvec;
          for jj=1:ly1
               for ii=1:lx1
                    pos = find(gno == El(elno).lglmap(ii,jj));
                    El(elno).soln(ii,jj) = soln(pos);
                    El(elno).un(ii,jj) = soln(pos);
                    posx = (jj-1)*lx1 + ii;
                    El(elno).unvec(posx,1) = soln(pos);
               end
          end
          if ifplot
               surf(El(elno).xm1,El(elno).ym1,El(elno).soln, 'EdgeColor', 'none');
               colorbar;
               hold on
%              view([0 90]);
          end
     end
     Elout = El;
     hold off
