function [xi_o,x0_o,vin_o,vout_o,vol1_o,vold_o,rnorm_o,ifinit_o,ik_o] = SpSnap(xi,x0,i,vol1,vold,vin,vout,ifsnap,sfreq,ifinit,ik,skryl,vlen)
% SPSNAP

   rnorm = norm(xi);

   if (~ifsnap) 
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
   end  

%  IFSNAP

   if i==1
     vol1  = zeros(vlen,1);
     vold  = zeros(vlen,1);
     vin(:,1)=zeros(vlen,1);
     vout(:,1)=zeros(vlen,1);
    
     vol1  = xi;
     x0    = xi;
     ifinit = 0;
   else
     if (~ifinit)
       vin(:,1)  = xi;
       vout(:,1) = zeros(vlen,1);
       ik = 1;
       vol1 = xi;

       x0 = xi;

       ifinit = 1;

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
     end
   end 


%  Apply SNAP
   if mod(i,sfreq)==0
     if (ifinit)

        dv = (xi - vol1);
        vout(:,ik) = dv;          % vout = (A-I)x

        if ik==1
%         Nothing to be done
          x0 = xi;
          vol1 = xi;
        else

          [Q,R]=qr(vout(:,1:ik),0);               % economy size

          [U,S,V] = svd(R);       % gives V nor VT
          U = Q*U;
          s = diag(S);
          [smin ind] = min(s);
          v1  = V(:,ind);
          us1 = U(:,ind)*smin;
          x1  = vin(:,1:ik)*v1;

          x0  = x1 + us1;
          vol1 = x0;

%         New residual
          rnorm = smin;
         
        end  

        ik = ik+1; 

        if mod(ik,skryl)==0
%         Restart                
          vin(:,1)  = x0;
          ik = 1;
        else  
          vin(:,ik)  = x0;
        end

     end        % if ~ifinit
   else
     vold = (xi-vold);
     rnorm = norm(vold);
     vold = xi;
     x0 = xi;
   end          % if mod(i,sfreq)

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

