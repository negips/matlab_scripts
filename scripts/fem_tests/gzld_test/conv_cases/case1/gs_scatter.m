function vl = gs_scatter(v,N,nels,periodic)
%% DSSUM  - as nek calls it
   
   if nels==1
      vl = v;
      return;
   end   

   for nel = 1:nels
      st=(nel-1)*N + 1;
      en=(nel)*N + 1;
      vl(:,nel) = v(st:en);
   end
         
   return  
%    End of gs_scatter
