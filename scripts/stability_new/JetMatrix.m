function [A,B] = JetMatrix(N,alpha,beta,Re)

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
    yinf = 2.0; % 10*blh
    Jac = yinf/2;

    D1 = (1/Jac)*D1;
    D2 = (1/Jac^2)*D2;
    D3 = (1/Jac^3)*D3;
    D4 = (1/Jac^4)*D4;

    zi=sqrt(-1);

    % mean velocity

    k2  = alpha^2 + beta^2;
    Nos = N+1;
    Nsq = N+1;
    vec = [0:N]';
    yphys = yinf/2*(D0(:,2) + 0);    %cos(pi*vec/N);         % Wall normal points
    y0  = min(yphys);                       % y0 
    theta = yphys./0.1;
    u   = sech(theta).^2;         % U velocity profile (physical)
%    du  = -2*(sech(theta).^2)*tanh(theta);     % y derivative (physical)
%    ddu = D2*D0*u;                       % 2nd derivative (physical)
    du  = D1*D0*u;
    ddu = D2*D0*u;

%    figure(500)  
%    plot(yphys,u)
%    hold on
%    plot(yphys,du, 'r')
%    plot(yphys,ddu, 'k')      
%    ylabel('U,dU,ddU')
%    xlabel('y')

% %    legend({'U', 'dU', 'ddU'}, 'Location', 'Best') 
%     legend({'U'}, 'Location', 'Best');  
%     hold off            

    % set up Orr-Sommerfeld matrix
    B11 = D2 - k2*D0;
    A11 = -(D4 - 2*k2*D2 + (k2*k2)*D0)/(zi*Re);             % Viscous term
    A11 = A11 + alpha*(u*ones(1,length(u))).*B11;           % alpha*U
%    A11 = A11 + (v*ones(1,length(v))).*((D3 -k2*D1)/zi);    % 1/i*V0*d/dy
    A11 = A11 - alpha*(ddu*ones(1,length(ddu))).*D0;        % alpha*d2U/dy2  
    er  = -200*zi;
    A11 = [er*[D0(1,:); D1(1,:)]; A11(3:Nos-2,:); ...       % This should be boundary conditions 
           er*[D1(Nos,:); D0(Nos,:)]];
    B11 = [D0(1,:); D1(1,:); B11(3:Nos-2,:); ... 
           D1(Nos,:); D0(Nos,:)];

%     Changing boundary conditions
%    A11 = [A11(1:Nos-2,:); ...       % This should be boundary conditions 
%           er*[D1(Nos,:); D0(Nos,:)]];
%    B11 = [B11(1:Nos-2,:); ... 
%           D1(Nos,:); D0(Nos,:)];



    % set up Squire matrix and (cross-term) coupling matrix
    A21 = beta*(du*ones(1,length(u))).*D0(1:Nos,:);         % Coupling term
    A22 = -(D2-k2*D0)/(zi*Re);                              % Viscous term
    A22 = A22 + alpha*(u*ones(1,length(u))).*D0;            % alpha*U      
%    A22 = A22 + (v*ones(1,length(v))).*D1/(zi);             % 1/i*V*d/dx

    B22 = D0;
    A22 = [er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];       % This should be boundary conditions
    A21 = [zeros(1,Nos); A21(2:Nsq-1,:); zeros(1,Nos)];
  
%    dbstop in SuctionMatrix at 66 
    % combine all the blocks 
    A = [A11 zeros(Nos,Nsq); A21 A22];
    B = [B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];

