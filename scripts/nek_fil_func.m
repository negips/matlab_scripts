%% Build basis change matrix
%% Taken directly from nek

clear
clc
close all

N=20;
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
          pht(kj) = Lj(k) - Lj(k-2);
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

% boydstonodal*G*inv(boydstonodal)

%%
% matrix for filter forcing

%G = zeros(N+1);
%G(N+1,N+1) =0.02;
%G(N,N) = 0.5;
%G(N,N) =1;
%G(N-1,N-1) =1;

%G = eye(N+1);

kut=6;
wght=1.0;
nx=N+1;
filter_func=zeros(N+1,1);
diag = zeros(N+1,N+1);
k0=nx-kut;
for k=k0+1:nx
  kk=k+nx*(k-1);
  amp=wght*(k-k0)*(k-k0)/(kut*kut);        % quadratic growth
  diag(k,k) = 1-amp;
  filter_func(k) = amp;
end

h1=figure;  
plot(filter_func, 'LineWidth', 2)
ylabel('$\mathcal{H}$', 'Interpreter', 'Latex', 'FontSize', 18)
xlabel('$N$', 'FontSize', 18)
%title('Filter Function shape')
xlim([1 N+1])
filename='filter_shape';

SaveFig(h1,filename,'./plots',1) 


%fil_mat = boydstonodal*G*inv(boydstonodal)

