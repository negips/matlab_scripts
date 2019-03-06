function [xi_o,x0_o,vin_o,vout_o,vol1_o,vold_o,rnorm_o,ifinit_o,ik_o] = SpSnap(xi,x0,i,vol1,vold,vin,vout,ifsnap,sfreq,ifinit,ik,skryl,vlen)
% SPSNAP

   rnorm = norm(xi);   

   if (ifsnap)
     if i==1
       vol1  = zeros(vlen,1);
       vold  = zeros(vlen,1);
        
       vol1  = xi;
       ifinit = 0;
     end
   end


%  Apply SNAP
   if (ifsnap)
     if mod(i,sfreq)==0
       if (~ifinit)

          vin(:,1)=xi;
          vout(:,1)=zeros(vlen,1);
          ik = 1;
          vol1  = xi;        % not in code. Probably initialization bug

          x0 = xi;

          ifinit = 1;
       else
          dv = (xi - vol1);
          vout(:,ik) = dv;          % vout = (A-I)x

          if ik==1
%           Nothing to be done
            x0 = xi;
          else

            [U,S,V] = svd(vout(:,1:ik),'econ');       % gives V nor VT
            s = diag(S);
            [smin ind] = min(s);
            v1  = V(:,ind);
            us1 = U(:,ind)*smin;
            x1  = vin(:,1:ik)*v1;
%            dbstop in SpSnap at 47

            x0  = x1 + us1;

%           New residual
            rnorm = smin;
           
          end  

          ik = ik+1;  

          if mod(ik,skryl)==0
%           Restart                
            vin(:,1)  = x0;
            ik = 1;
          else  
            vin(:,ik)  = x0;
          end

%          dbstop in BoostConv at 46
         
%          [norm(xi) norm(x0) smin norm(vol1)]
          vol1 = x0;
          vold = x0;            % also not in code

%          dbstop in SpSnap at 72

       end
     else
       vold = (xi-vold);
       rnorm = norm(vold);
       vold = xi;
       x0 = xi;

     end          % if mod(i,sfreq)
   end            % ifsnap    


   xi_o=xi;
   x0_o=x0;
   vin_o=vin;
   vout_o=vout;
   vol1_o=vol1;
   vold_o=vold;
   rnorm_o=rnorm;
   ifinit_o=ifinit;
   ik_o=ik;


return

