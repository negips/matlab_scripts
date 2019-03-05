function [xi_o,x0_o,dv_o,vin_o,vout_o,vol1_o,vold_o,rnorm_o,ifinit_o,ib_o] = BoostConv(xi,x0,i,vol1,vold,vin,vout,ifboost,bfreq,ifinit,ib,bkryl,vlen)
% BoostConv

   dv = [];
   rnorm = 0.;
   c   = zeros(bkryl,1);

   if (ifboost)
     if i==1
       vol1 = zeros(vlen,1);
       vold  = zeros(vlen,1);
        
       vold  = xi;
       ifinit = 0;
     end
   end


%  Apply Boostconv
   if (ifboost)
     if mod(i,bfreq)==0
       if (~ifinit)
          ib = 1;
          dv = (xi - vol1);

          vin(:,1)=dv;
          vout(:,1)=dv;
          ifinit = 1;
          
          vol1  = xi;        % not in code. Probably initialization bug
       else
          
          dv = (xi - vol1);

          vout(:,ib) = vout(:,ib)-dv;             % vout = rn_1 - rn
          vin(:,ib)  = vin(:,ib)-vout(:,ib);      % vin  = rn

          c_bc  = vout(:,1:ib)'*dv;
          dd_bc = vout(:,1:ib)'*vout(:,1:ib);
          c2    = dd_bc\c_bc;               % invert dd_bc
          c(1:length(c2),1) = c2;

          ib = mod(ib,bkryl)+1;
          vout(:,ib) = dv;                  % vout+1 = rn

%          dbstop in BoostConv at 46

%         New residual
          dv = dv + vin*c;
          vin(:,ib) = dv;

          rnorm = norm(dv);
         
          x0 = vol1+dv;
          vol1 = x0;
          vold = x0;            % also not in code

       end
     else
       vold = (xi-vold);
       rnorm = norm(vold);
       vold = xi;

     end          % if mod(i,bfreq)
   end            % ifboost     

   xi_o=xi;
   x0_o=x0;
   dv_o=dv;
   vin_o=vin;
   vout_o=vout;
   vol1_o=vol1;
   vold_o=vold;
   rnorm_o=rnorm;
   ifinit_o=ifinit;
   ib_o=ib;


return

