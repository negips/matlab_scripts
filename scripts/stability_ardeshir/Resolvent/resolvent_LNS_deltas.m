% This script runs the temporal 1D stability problem
%
% (c) Ardeshir Hanifi 2014
%
clear all
%% Numerical parameters
N = 300; % number of Cehebychev modes

%% Flow parameters
Re=400;
alpha=0;
%beta=0.45;
omegas=0;
beta=linspace(1e-2,1,100);
% beta=2; %(0.01:0.01:0.3);
% omegas=0.0;
flowcase=2; % 1: channel flow, 2: BL

%% Select the collocation nodes and set up derivative operator
[y, D] = chebdif(N, 2);

%% Velocity profile corresponding to Plane Poiseuille Flow
if (flowcase==1)
    U.u=1-y.^2;
    U.uy=-2*y;
    U.w=zeros(size(y)); % This field is required as the code is written for a 3D flow
    U.wy=zeros(size(y)); % This field is required as the code is written for a 3D flow
    W = iwt(N);
    ybl=y;
else
    sweep = 0.;
    deltas=1.7212;
    mfac=0;
    ymax=30;
    ybl=(y+1)*ymax/2*deltas;
    U=similar_fsc(ybl,mfac);
    for k=1:2
        D(:,:,k)=D(:,:,k)*(2/ymax)^k;
    end
    U.uy =D(:,:,1)*U.u;
    U.uyy=D(:,:,2)*U.u;
    W = iwt(N)*ymax/2;
    ybl=ybl/deltas;
    we=tand(sweep);
    U.w = U.w*we;
    U.wy =D(:,:,1)*U.w;
    U.wyy=D(:,:,2)*U.w;
 end
%% Set up the eigenvalue problem

G=zeros(length(omegas),length(beta));

for k2=1:length(beta)
    
[A,B] = LNS_temp(U, D, alpha, beta(k2), Re);

%% Solve the eigenvalue problem
[eigvec, eigval] = eig(A,B);
eigval=diag(eigval);
TF=isfinite(eigval); indeg=(TF~=0);
eigval=eigval(indeg);
eigvec=eigvec(:,indeg);
nmode=length(eigval);
%% Plot spectrum and eigenvectors
% header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
%               {', beta= '},{num2str(beta)});
% plot_LNS(eigval,eigvec,y,header)
% figure(20),clf
% plot(real(eigval(1:nmode)),imag(eigval(1:nmode)),'o')
% ylim([-0.1,0.05]); xlim([0,(sqrt(alpha^2+beta^2))])
% grid on

%% 

[eigvec, M]=setM( eigvec, W );

F=sqrtm(M);

%% Resolvant norm

for k=1:length(omegas)
    FwF=F*diag(1./(eigval-omegas(k)))/F;
    G(k,k2)=max(abs(svd(FwF)));
    disp(['Omega: ',num2str(omegas(k)),' Beta: ',num2str(beta(k2)),' G: ',num2str(G(k,k2))])
end

end 
%%
figure(10); clf
semilogy(omegas,G(:,1),'o-r')
xlabel('omega'); ylabel('G')
str=['Re= ',num2str(Re),', alpha= ',num2str(alpha),', beta= ',num2str(beta(1))];
title(str)

figure(11); clf
plot(beta,G(1,:),'o-r')
xlabel('beta'); ylabel('G')
str=['Re= ',num2str(Re),', alpha= ',num2str(alpha),', omega= ',num2str(omegas(1))];
title(str)
%%
FwF=F*diag(1./(eigval-omegas(1)))/F;
[u1,s,v1]=svd(FwF);
%[v1,eg]=eig(B'*B); eg=diag(eg); [~,indeg]=sort(eg,'descend'); v1=v1(:,indeg); eg=eg(indeg);
K0=F\v1(:,1);
Kt=F\u1(:,1)*s(1,1);
vopt=eigvec*K0;
vmax=eigvec*Kt;

figure(31);plot(abs(vopt(1:N)),ybl,abs(vopt(N+1:2*N)),ybl,...
        abs(vopt(2*N+1:3*N)),ybl);
legend('u','v','w');title('optimal initia')
eint=vopt'*blkdiag(W,W,W,W*0)*vopt;
disp(['Initial energy: ',num2str(abs(eint))]);

figure(32);plot(abs(vmax(1:N)),ybl,abs(vmax(N+1:2*N)),ybl,...
        abs(vmax(2*N+1:3*N)),ybl);
legend('u','v','w'); title('optimal final')
emax=vmax'*blkdiag(W,W,W,W*0)*vmax;
disp(['Optimal energy: ',num2str(abs(emax))]);