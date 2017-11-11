% Transient Growth
%
% compute the Orr-Sommerfeld matrix for three-
% dimensional Poiseuille or Couette flows,  
% compute the energy weight matrix and determines 
% transient growth curve G(t)
%
% INPUT 
%
% Re        = Reynolds number
% alpha     = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)  
% N         = total number of modes for normal velocity
% T         = compute maximum growth in time interval [0 T]
%

%    clear 
%    clc
    close all    
    
    global D0 D1 D2 D4 
    global qb
    global yphys  

    zi = sqrt(-1);

    %...input data
%    iflow  = input('Poiseuille (1) or Couette flow (2) or Asymptotic Suction (3)');
%    if (iflow == 3)
%      V0 = input('Suction velocity');
%    end    
%    N      = input('Enter the number of Chebyshev polynomials: ');
%    Re     = input('Enter the Reynolds number: ');
%    alpha  = input('Enter alpha: ');
%    beta   = input('Enter beta: ');
%    Tmax   = input('Enter Tmax: '); 
    iflow=3;
    N=1024;
%    Re=54385;
    Re=500;
    alpha=1.0;
    beta=1.0;
    Tmax=100;

    T      = [0 Tmax];
    %...generate Chebyshev differentiation matrices
    [D0,D1,D2,D3,D4] = ChebMat2(N);

    %...set up Orr-Sommerfeld matrices A and B 
    if (iflow == 1)
      [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
    elseif (iflow == 2)
      [A,B] = CouetteMatrix(N,alpha,beta,Re);
    else
      [A,B] = SuctionMatrix(N,alpha,beta,Re);
    end

    %...generate energy weight matrix
    k2 = alpha^2 + beta^2;
    M  = EnergyMatrix(N+1,N+1,k2);

    %...compute the Orr-Sommerfeld matrix (by inverting B)
    OS = inv(B)*A;

    % compute and plot the eigenvalue spectrum
    eee = eig(OS); 
    figure(4); 
    plot(real(eee),imag(eee),'o')
    axis([-.1 alpha*1.1 -1 .1]); 
    grid on
    
    %...compute the optimal
    [flowin,flowot,gg] = Optimal(OS,T,M,k2,1);

    %...graphics
    figure(2) 
    semilogy(gg(:,1),gg(:,2));
    grid on
%     figure(3)   
%     plot(gg(:,1),gg(:,2));
%     grid on
