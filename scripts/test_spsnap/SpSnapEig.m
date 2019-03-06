function [xi_o,x0_o,vin_o,vout_o,wout_o,vol1_o,vold_o,rnorm_o,ifinit_o,ik_o] = SpSnapEig(xi,x0,i,vol1,vold,vin,vout,wout,ifsnap,sfreq,ifinit,ik,skryl,vlen)
% SPSNAP

   rnorm = norm(xi);

   if (~ifsnap) 
     xi_o=xi;
     x0_o=x0;
     vin_o=vin;
     vout_o=vout;
     wout_o=wout;
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
       ik = 1;
       vol1 = xi;

       x0 = xi;

       ifinit = 1;

       xi_o=xi;
       x0_o=x0;
       vin_o=vin;
       vout_o=vout;
       wout_o=wout;
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

        beta = vol1'*vol1;
        vortho = vol1*(vol1'*xi)/beta;
        vin(:,ik)  = xi;              % vin  = Ax
        vout(:,ik) = vortho;          % vout = (xx^T/||x||^2)Ax
        wout(:,ik) = xi - vortho;     % wout = (I - xx^T/||x||^2)Ax

%        dbstop in SpSnapOrtho at 66

        if ik==1
%         Nothing to be done
          x0 = xi;
          vol1 = xi;
        else

          [U,S,V] = svd(wout(:,1:ik),'econ');       % gives V nor VT
          s = diag(S);
          [smin ind] = min(s);
          v1  = V(:,ind);
          us1 = U(:,ind)*smin;
          x1  = vin(:,1:ik)*v1;

%          dbstop in SpSnapOrtho at 79          

          x0  = x1;           % New start vector
          vol1 = x0;          

%         New residual
          rnorm = smin;
         
        end  

        ik = ik+1; 

        if mod(ik,skryl)==0
%         Restart                
          ik = 1;
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
   wout_o=wout;
   vol1_o=vol1;
   vold_o=vold;
   rnorm_o=rnorm;
   ifinit_o=ifinit;
   ik_o=ik;


return

