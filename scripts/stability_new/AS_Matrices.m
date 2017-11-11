function [A,B] = AS_Matrices(N,alpha,beta,Re,Uinf,V0)
%
% la*B*q=A*q
%
% N        = number of modes
% alpha    = alpha
% beta     = beta 
% Re       = Reynolds number
% Uinf     = Free stream velocity
% V0       = Suction velocity
%
% Erik Bostr?m
% 20170101
%
    global D0 D1 D2 D4
    
    zi=1i;

    k2  = alpha^2 + beta^2;
    Nos = N+1;
    Nsq = N+1;
    vec = (0:N)';
    one = ones(1,length(vec));
    y   = (cos(pi*vec/N)+1)/2;
    u   = Uinf*(one'-exp(-abs(V0)*Re.*y));
    du  = abs(V0)*Re*exp(-abs(V0)*Re.*y);
    d2u = -(abs(V0)*Re)^2*exp(-abs(V0)*Re.*y);
    O   = zeros(Nos,Nsq);
    
    B11 = D2 - k2*D0;
    B22 = D0;
    
    A11 = (alpha*u*one + V0*D1/zi).*B11 ...
          + alpha*d2u*one.*D0 ...
          - (D4 -2*k2*D2 + (k2*k2)*D0)/(zi*Re);
        
    er  = -200*zi; % Why -200??

    A11 = [ er*D0(1,:)      ; 
            er*D1(1,:)      ; 
            A11(3:Nos-2,:)  ; 
            er*D1(Nos,:)    ; 
            er*D0(Nos,:)    ];
    
    A21 = beta*du*one.*D0(1:Nos,:);
    A21 = [ zeros(1,Nos)    ; 
            A21(2:Nsq-1,:)  ; 
            zeros(1,Nos)]   ;
    
    A22 = alpha*u*one.*D0 + V0*D1/zi - (D2-k2*D0)/(zi*Re);
    A22 = [ er*D0(1,:)      ; 
            A22(2:Nsq-1,:)  ; 
            er*D0(Nsq,:)    ];

        
    A = [ A11  O   ; 
          A21  A22 ];

    B = [ B11  O   ; 
          O    B22 ];

