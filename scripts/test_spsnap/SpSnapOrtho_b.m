function [xi_o,x0_o,X_o,Y_o,Z_o,W_o,vol1_o,vold_o,rnorm_o,ifinit_o,ik_o] = SpSnapOrtho(xi,x0,b,i,vol1,vold,X,Y,Z,W,ifsnap,sfreq,ifinit,ik,skryl,vlen)
% SPSNAP

   ifthick = 1;
   nthick = 5; 

   rnorm = norm(xi);

   if (~ifsnap) 
     xi_o=xi;
     x0_o=x0;
     X_o=X;
     Y_o=Y;
     Z_o=Z;
     W_o=W;
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
     Y(:,1)=zeros(vlen,1);
     Z(:,1)=zeros(vlen,1);
    
     vol1  = xi;
     ifinit = 0;
   else
     if (~ifinit)
       ik = 1;
       vol1 = xi;

       x0 = xi;

       ifinit = 1;

       xi_o=xi;
       x0_o=x0;
       X_o=X;
       Y_o=Y;
       Z_o=Z;
       W_o=W;
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

        beta = (b'*b);

        y = xi;
        z = b*(b'*xi)/beta;
        w = y - z;
        Y(:,ik) = y;         % Y  = Ax
        Z(:,ik) = z;         % Z = (bb^T/||b||^2)Ax
        W(:,ik) = y - z;     % W = (I - xx^T/||x||^2)Ax = EAx

        if ik==1
%         Nothing to be done
          x0 = xi;
          vol1 = xi;
        else

          [U,S,V] = svd(W(:,1:ik),'econ');       % gives V nor VT
          s = diag(S);
          [smin ind] = min(s);
          v1  = V(:,ind);
          us1 = U(:,ind)*smin;
          x1  = Y(:,1:ik)*v1;

          x0  = (b'*b)/(b'*x1)*x1;
%          x0 = x1;  
          vol1 = x0;          

%         New residual
          rnorm = smin;
         
        end  

        ik = ik+1; 

        if ik>skryl
%         Restart                
          ik = 1;
          if (ifthick)

            vthick = V(:,end-nthick+1:end);
%            dbstop in SpSnapOrtho_b at 103

            Y(:,1:nthick) = Y*vthick;
            W(:,1:nthick) = W*vthick;
            ik=nthick+1;
          end  
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
   X_o=X;
   Y_o=Y;
   Z_o=Z;
   W_o=W;
   vol1_o=vol1;
   vold_o=vold;
   rnorm_o=rnorm;
   ifinit_o=ifinit;
   ik_o=ik;


return

