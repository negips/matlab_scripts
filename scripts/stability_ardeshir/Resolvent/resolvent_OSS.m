% This script runs the temporal 1D stability code
%
% (c) Ardeshir Hanifi 2014
%
clear all
%% Numerical parameters
N = 300; % number of Cehebychev modes

%% Flow parameters
Re=300;
alpha=0;
beta=0.45;
omegas=0;
%beta=(0.1221);
%omegas=(-0.06:0.001:0.06);
%beta=(0.01:0.1:2);
%omegas=0.0;

flowcase=2; % 1: channel flow, 2: BL

%% Select the collocation nodes and set up derivative operator
[y, D] = chebdif(N, 4);

%% Velocity profile corresponding to Plane Poiseuille Flow
if (flowcase==1)
    U.u=1-y.^2;
    U.uy=-2*y;
    U.uyy=-2+y*0;
    W = iwt(N);
else
    mfac=0;
    ymax=100;
    ybl=(y+1)*ymax/2;
    U=similar_fsc(ybl,mfac);
    for k=1:4
        D(:,:,k)=D(:,:,k)*(2/ymax)^k;
    end
    U.uy =D(:,:,1)*U.u;
    U.uyy=D(:,:,2)*U.u;
    W = iwt(N)*ymax/2;
end

%% Set up the eigenvalue problem

G=zeros(length(omegas),length(beta));

for k2=1:length(beta)
[A,B] = OSS_temp(U, D, alpha, beta(k2), Re);

%% Solve the eigenvalue problem
[eigvec, eigval] = eig(A,B); eigval=diag(eigval);

% TF=isfinite(eigval); indeg=(TF~=0);
[~,indeg]=sort(-imag(eigval));
eigval=eigval(indeg);
eigvec=eigvec(:,indeg);
nmode=length(eigval)-10;
eigval=eigval(1:nmode);
eigvec=eigvec(:,1:nmode);

%%
[ eigvec, M ]=setMoss( eigvec, W, D, alpha, beta(k2));

F=sqrtm(M);

%% Resolvant norm

for k=1:length(omegas)
    FwF=F*diag(1./(eigval-omegas(k)))/F;
    G(k,k2)=max(abs(svd(FwF)));
    disp(['Omega: ',num2str(omegas(k)),' Beta: ',num2str(beta(k2)),' G: ',num2str(G(k,k2))])
end

end 
%%
figure(20); clf
semilogy(omegas,G(:,1),'o-r')
xlabel('omega'); ylabel('G')
str=['Re= ',num2str(Re),', alpha= ',num2str(alpha),', beta= ',num2str(beta(1))];
title(str)

figure(21); clf
plot(beta,G(1,:),'o-r')
xlabel('beta'); ylabel('G')
str=['Re= ',num2str(Re),', alpha= ',num2str(alpha),', omega= ',num2str(omegas(1))];
title(str)
%%
FwF=F*diag(1./(eigval-omegas(1)))/F;
[u1,s,v1]=svd(FwF);
K0=F\v1(:,1);
Kt=F\u1(:,1)*s(1,1);
vopt=eigvec*K0;
vmax=eigvec*Kt;

k2=alpha*alpha+beta*beta;
Dv=D(:,:,1)*vopt(1:N);
u=(alpha*Dv-beta*vopt(N+1:2*N))/(-1i*k2);
w=(beta*Dv+alpha*vopt(N+1:2*N))/(-1i*k2);
figure(41);plot(abs(u),ybl,'r',abs(vopt(1:N)),ybl,'b',abs(w),ybl,'k');
legend('u','v','w');title('optimal initia')

Dv=D(:,:,1)*vmax(1:N);
u=(alpha*Dv-beta*vmax(N+1:2*N))/(-1i*k2);
w=(beta*Dv+alpha*vmax(N+1:2*N))/(-1i*k2);
figure(42);plot(abs(u),ybl,'r',abs(vmax(1:N)),ybl,'b',abs(w),ybl,'k');
legend('u','v','w'); title('optimal final')
