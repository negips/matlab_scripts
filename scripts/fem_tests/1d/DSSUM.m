function var = DSSUM(v,nels,periodic)
%% DSSUM  - as nek calls it
   
   var = v;    
   if nels==1
     if periodic
          v(1) = (v(1) + v(end));
          v(end) = v(1);
          var = v;
          return;
     else
          return;
     end 
   end   

   for nel = 1:nels
     if nel==1
          v(end,nel) = v(end,nel)+v(1,nel+1);

          if periodic
               v(1,1) = (v(1,1) + v(end,nels));
          end
     elseif nel<nels
          v(1,nel) = v(end,nel-1);
          v(end,nel) = v(end,nel) + v(1,nel+1);

     elseif nel==nels
          v(1,nel) = v(end,nel-1);

          if periodic
               v(end,nel) = v(1,1);
          end
     end
   end         
   var = v;  
   return  
%    End of DSSUM
