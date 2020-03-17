function [spectonodal nodaltospec] = Leg2Nodal(N)


[x w] = lglnodes(N);
x = x(end:-1:1);
nx = N+1;
kj = 0;
n = nx-1;
for j = 1:nx
     z=x(j);
     Lj = legendrePoly(n,z);
     kj = kj+1;
     pht(kj) = Lj(1);
     kj = kj+1;
     pht(kj) = Lj(2);
     for k = 3:nx
          kj = kj+1;
          pht(kj) = Lj(k);
     end
end

pht = reshape(pht,nx,nx);
spectonodal = transpose(pht);
phi = pht;

nodaltospec = inv(spectonodal);

