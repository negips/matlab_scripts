function Elout = vmult(El,nelv,lx1,ly1)

     for ee =1:nelv
          vmult = 0*El(ee).lglmap; 

          for jj=1:ly1
          for ii=1:lx1
               glno = El(ee).lglmap(ii,jj);
               reps = 0;
               for e2 =1:nelv
                    ind = El(e2).lglmap == glno;
                    reps = reps+sum(ind(:));
               end
               vmult(ii,jj) = reps;
          end
          end
          El(ee).vmult = vmult;
     end

     Elout = El;
     return
          
