function [L] =legendrePoly(N,x)

if N==0
     L(1) = 1;
     return
end

if N==1
     L(1) = 1;
     L(2) = x;
     return
end

L = zeros(N+1,1);
L(1) = 1;
L(2) = x;

for i = 3:N+1
     j = i-1;
     L(i) = ((2*j-1)*x*L(i-1) - (j-1)*L(i-2))/j;
end

end
