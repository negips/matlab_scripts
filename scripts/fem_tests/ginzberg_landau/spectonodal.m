
%% Build basis change matrix
%% Taken directly from nek
addpath '../'

N = 5;
x = lglnodes(N);
x=x(end:-1:1);

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
          pht(kj) = Lj(k);% - Lj(k-2);
     end
end

pht = reshape(pht,nx,nx);
spec2nodal = transpose(pht);
phi = pht;

boyd = 0;
if boyd
     for i = 3:nx
     for j = 1:nx
          phi(i,j) = pht(i,j) - pht(i-2,j);
     end
     end
end
boyds2nodal = transpose(phi);
