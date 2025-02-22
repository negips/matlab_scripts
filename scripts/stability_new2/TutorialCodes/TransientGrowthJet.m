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

    clear 
    
    global D0 D1 D2 D4 
    global qb
    
    zi = sqrt(-1);

    %...input data
%    iflow  = input('Poiseuille (1) or Couette flow (2) ');
%    N      = input('Enter the number of Chebyshev polynomials: ');
%    Re     = input('Enter the Reynolds number: ');
%    alpha  = input('Enter alpha: ');
%    beta   = input('Enter beta: ');
%    Tmax   = input('Enter Tmax: ');

    iflow=3;
    N=200;
    Re=10.0;
    alpha=0.720;
    beta =0.0;
    Tmax =50;

    T      = [0 Tmax];
    %...generate Chebyshev differentiation matrices
    [D0,D1,D2,D4] = ChebMat(N);

    %...set up Orr-Sommerfeld matrices A and B 
    if (iflow == 1)
      [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
    elseif (iflow==2)
      [A,B] = CouetteMatrix(N,alpha,beta,Re);
    elseif (iflow==3)
      [A,B] = JetMatrix2(N,alpha,beta,Re);
    end

    %...generate energy weight matrix
    k2 = alpha^2 + beta^2;
    M  = EnergyMatrix(N+1,N+1,k2);

    %...compute the Orr-Sommerfeld matrix (by inverting B)
    OS = inv(B)*A;  
    % compute and plot the eigenvalue spectrum
    eee = eig(OS); 
    figure(3); 
    plot(real(eee),imag(eee),'o')
    axis([-.1 alpha*1.1 -1 .1]); 
    grid on
    
    %...compute the optimal
    [flowin,flowot,gg] = Optimal(OS,T,M,k2,1);

    %...graphics
    figure(1) 
    semilogy(gg(:,1),gg(:,2));
    grid on
    figure(2)   
    plot(gg(:,1),gg(:,2));
    grid on
