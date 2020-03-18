function [A,B] = JetMatrix2(N,alpha,beta,Re);

%
% Function to create Orr-Sommerfeld matrices using Chebyshev 
% pseudospectral discretization for plane Couette flow 
% profile
%
% N   = number of even or odd modes
% alpha     = alpha
% Re       = Reynolds number

    global D0 D1 D2 D3 D4
    global yphys  

    blh = -log(0.01);
    yinf = 5;
    Jac = yinf;

    D1 = (1/Jac)*D1;
    D2 = (1/Jac^2)*D2;
    D3 = (1/Jac^3)*D3;
    D4 = (1/Jac^4)*D4;

    zi=sqrt(-1);

%   mean velocity

    k2  = alpha^2 + beta^2;
    Nos = N+1;
    Nsq = N+1;
    vec = [0:N]';

    yphys = yinf*D0(:,2);       %cos(pi*vec/N);         % Wall normal points
    theta = yphys; 

    u   = sech(theta).^2;         % U velocity profile (physical)
%    du  = D1*D0*u;
%    ddu = D2*D0*u;
    du  = -2*tanh(theta).*sech(theta).*sech(theta);
    ddu =  4*(tanh(theta).*sech(theta)).^2 - 2*sech(theta).^4;

    figure(500)
    plot(yphys,ddu)

    % set up Orr-Sommerfeld matrix
    B11 = D2 - k2*D0;
    A11 = -(D4 - 2*k2*D2 + (k2*k2)*D0)/(zi*Re);
    A11 = A11 + alpha*(u*ones(1,length(u))).*B11;
    A11 = A11 - alpha*(ddu*ones(1,length(ddu))).*D0;
    er  = -2000*zi;
    A11 = [er*[D0(1,:); D1(1,:)]; A11(3:Nos-2,:); ... 
           er*[D1(Nos,:); D0(Nos,:)]];
    B11 = [D0(1,:); D1(1,:); B11(3:Nos-2,:); ... 
           D1(Nos,:); D0(Nos,:)];

    % set up Squire matrix and (cross-term) coupling matrix
    A21 = beta*(du*ones(1,length(u))).*D0;
    A22 = alpha*(u*ones(1,length(u))).*D0-(D2-k2*D0)/(zi*Re);
    B22 = D0;
    A22 = [er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
    
    % combine all the blocks 
    A = [A11 zeros(Nos,Nsq); A21 A22];
    B = [B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];

