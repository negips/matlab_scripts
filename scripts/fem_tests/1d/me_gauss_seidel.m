function [soln err err_decay itrs] = me_gauss_seidel(A,x0,b,maxit,tol)


[nx ny nz] = size(A);

L = 0*A;
U = L;

nels = nz;

x0 = zeros(nx,nels);
loc_x = zeros(nx,nels);
err_decay = zeros(maxit,1);

for nel=1:nels
     [L(:,:,nel) U(:,:,nel)] = lu(A(:,:,nel));             % LU decomposition
end

% xk = (L^-1)(b - ux);
conv = 0;
itrs =0;
while (~conv) && (itrs<maxit)
itrs=itrs+1;
for nel=1:nels
     loc_x(:,nel) = inv(L(:,:,nel))*(b(:,nel) - U(:,:,nel)*loc_x(:,nel));
     if nel~=nels
          loc_x(1,nel+1) = loc_x(end,nel);
     else
          loc_x(1,1) = loc_x(1,1);
     end
end

for nel=1:nels
     r(:,nel) = b(:,nel) - A(:,:,nel)*loc_x(:,nel);
     err(nel) = norm(r);

     if norm(b(:,nel))~=0
          rerr(nel)=err(nel)/norm(b(:,nel));
     else
          rerr(nel)=err(nel);
     end

     ifconv(nel) = 0;
     if rerr(nel)<tol
          ifconv(nel) = 1;
     end
end

if sum(ifconv)==nels
     conv=1;
end

err_decay(itrs) = sum(rerr);
end

soln=loc_x;
err=rerr; 
return



