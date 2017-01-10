%    Testing some numbers with Boyd mode transformation

clear
clc
close all

N = 8;         % polynomial order

Un = rand(N+1,1);

if (N<2)
     display('N must be greater than 2!');
     break;
end


boyd_tr = zeros(N+1);
boyd_tr(1,1) = 1;
boyd_tr(2,2) = 1;


for i = 3:N+1
     boyd_tr(i,i) = 1;
     boyd_tr(i,i-2) = -1; 
end
boyd_inv = inv(boyd_tr);

G = eye(N+1);

wght = [0.0];      % In reverse order wght(N+1,N,N-1,....1)

l1 = length(wght);

for i = 0:l1-1
     ind = N+1-i;
     G(ind,ind) = wght(i+1);
end

Utr = boyd_inv*G*boyd_tr;

boyd_tr;
boyd_inv;
G;
Utr;


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
boyd = 1;

if boyd
     for i = 3:nx
     for j = 1:nx
          phi(i,j) = pht(i,j) - pht(i-2,j);
     end
     end
end

boydstonodal = transpose(phi);


for j = 1:100000

Un = 2*rand(N+1,1)-1;
Un = 0.5*(2*rand(N+1,1)-1 + Un(end:-1:1));

legspec = inv(spectonodal)*Un;
nodalvals = boydstonodal*G*inv(boydstonodal)*Un;
legspec2 = inv(spectonodal)*nodalvals;

modorg(j) = legspec(N-1);
modnew(j) = legspec2(N-1);
deltae(j) = modnew(j)^2 - modorg(j)^2;
modnorg(j) = legspec(N+1);
tdeltae(j) = deltae(j) + legspec2(N+1)^2 - legspec(N+1)^2;
end

%disp('Spectra');
%display([legspec legspec2]);
%display('Nodal vals');
%display([Un nodalvals]);

scatter(modorg./modnorg,tdeltae)
xlim([-10 10])

figure
scatter(modnorg,modorg, 'r')



