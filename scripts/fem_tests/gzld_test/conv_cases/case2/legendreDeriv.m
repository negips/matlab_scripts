function [dL] =legendreDeriv(N,x)

if N==0
  dL(1) = 0;
  return
end

if N==1
  dL(1) = 0;
  dL(2) = 1;
  return
end

dL = zeros(N+1,1);
dL(1) = 0;
dL(2) = 1;

L = legendrePoly(N,x);

for i = 3:N+1
  n = i-1;
  dL(i) = n*L(n) + x*dL(n);
end

return
