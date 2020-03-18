function Elout = BdryIndex(El,nelv,lx1,ly1);

     for elno = 1:nelv
          cnt=0;
          bdryind = zeros(2*lx1 + 2*ly1 -4,2);
          bdrygno = zeros(2*lx1 + 2*ly1 -4,1);
          for jj=1:ly1
          for ii=1:lx1
               if El(elno).vmult(ii,jj)>1      % Connected node
                    cnt=cnt+1;
                    bdryind(cnt,1) = ii;
                    bdryind(cnt,2) = jj;
                    bdrygno(cnt) = El(elno).lglmap(ii,jj);
               end
          end
          end
          El(elno).bdryind = bdryind;
          El(elno).bdrygno = bdrygno;
     end
     Elout = El;
     return
