function [El ifconv] = SolveSem2D(El,lx1,ly1,nelv,max_it,tol)

%  -- Iterative solve routine --

% for Multi element SEMs.     

% build and subtract forcing on boundary nodes

bdry_len = 2*lx1+2*ly1-4;
ifconv=0;
h=figure;
for iter=1:max_it
     
     for elno=1:nelv
          El(elno).scrtch1 = El(elno).b;
          for e2 =1:nelv
               if e2==elno
                    continue;
               end

               for kk=1:bdry_len
                    glno = El(elno).bdrygno(kk);
                    ind = find(El(e2).bdrygno == glno);
                    if ~isempty(ind)
                         ix = El(e2).bdryind(ind,1);
                         iy = El(e2).bdryind(ind,2);
                         b2 = El(e2).A(ix,:,iy)*El(e2).soln(:,iy);

                         ix = El(elno).bdryind(kk,1);
                         iy = El(elno).bdryind(kk,2);
                         El(elno).scrtch1(ix,iy) = El(elno).b(ix,iy) - b2; 
                    end
               end
          end

%         Either build all bdry forcing for elements and then solve.
%         Or solve for each element as soon as we build the forcing vector.
          for jj=1:ly1
               [El(elno).soln(:,jj) cflag iter relres] = pcg(El(elno).A(:,:,jj),El(elno).scrtch1(:,jj));
          end
     end

%     solnplot(El,nelv,h)
    
     El = DSAV2D_soln(El,lx1,ly1,nelv);

     ifconv = CheckConv(El,lx1,ly1,nelv,tol);
     if ifconv==1
          break
     end

end

close(h);

return

%-------------------------------------------------- 
function El = DSAV2D_soln(El,lx1,ly1,nelv)
     
bdry_len = 2*lx1+2*ly1-4;
     
for elno=1:nelv
     El(elno).scrtch2 = El(elno).soln;
     for e2 =1:nelv
          if e2==elno
               continue;
          end

          for kk=1:bdry_len
               glno = El(elno).bdrygno(kk);
               ind = find(El(e2).bdrygno == glno);
               if ~isempty(ind)
                    ix = El(e2).bdryind(ind,1);
                    iy = El(e2).bdryind(ind,2);
                    b2 = El(e2).soln(ix,iy);
                    
                    ix = El(elno).bdryind(kk,1);
                    iy = El(elno).bdryind(kk,2);
                    El(elno).scrtch2(ix,iy) = El(elno).scrtch2(ix,iy) + b2; 
               end
          end
     end
end 

for elno=1:nelv
     El(elno).soln = El(elno).scrtch2./El(elno).vmult;
end

return
%-------------------------------------------------- 

function ifconv = CheckConv(El,lx1,ly1,nelv,tol)
     
     lr_e=zeros(nelv,1);
     for elno=1:nelv
          r = zeros(lx1,ly1);
          for jj=1:ly1
               r(:,jj) = El(elno).A(:,:,jj)*El(elno).soln(:,jj) - El(elno).bl(:,jj);
          end
          r = r.^2;
          lr_e(elno) = sqrt(sum(sum(r)));
          if lr_e(elno)>tol
               ifconv=0;
               return
          end
     end
     ifconv=1;
     return

%-------------------------------------------------- 







