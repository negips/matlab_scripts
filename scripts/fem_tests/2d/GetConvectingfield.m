function [Cx Cy] = GetConvectingfield(xm1,ym1)


Cx=zeros(size(xm1));
Cy=zeros(size(ym1));

n = 1;

alpha=n*2*pi;
beta=n*2*pi;

m=10;
alpha2=m*2*pi;
beta2=m*2*pi;


V0=0.5;
U0=1.0;

du=0.1*U0;
dv=-du*alpha/beta;

Cx = U0 + du*sin(alpha*xm1 + beta*ym1) + 0.0*du*sin(beta2*ym1);
Cy = V0 + dv*sin(alpha*xm1 + beta*ym1) + 0.0*dv*sin(alpha2*xm1);

% Might need to project lower polynomial order.
% Does it still remain divergence free?

