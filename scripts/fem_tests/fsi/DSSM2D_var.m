function El = DSSM2D_var(El,nelv,varname)


for elno=1:nelv
  evalstr = ['El(' num2str(elno) ').scrtch1 = El(' num2str(elno) ').' varname ';'];
  eval(evalstr)
  evalstr = ['El(' num2str(elno) ').scrtch2 = El(' num2str(elno) ').' varname ';'];
  eval(evalstr)
end


[lx1 ly1] = size(El(1).scrtch1);
     
bdry_len = 2*lx1+2*ly1-4;
     
for elno=1:nelv
%     El(elno).scrtch1 = El(elno).bl;
     for e2 =1:nelv
          if  e2==elno
               continue;
          end

          for kk=1:bdry_len
               glno = El(elno).bdrygno(kk);
               ind = find(El(e2).bdrygno == glno);
               if ~isempty(ind)
                    ix = El(e2).bdryind(ind,1);
                    iy = El(e2).bdryind(ind,2);
                    b2 = El(e2).scrtch2(ix,iy);
                    
                    ix = El(elno).bdryind(kk,1);
                    iy = El(elno).bdryind(kk,2);
                    El(elno).scrtch1(ix,iy) = El(elno).scrtch1(ix,iy) + b2; 
               end
          end
     end
end 

for elno=1:nelv
  evalstr=['El(' num2str(elno) ').' varname '= El(' num2str(elno) ').scrtch1;'];
  eval(evalstr)
end


return

