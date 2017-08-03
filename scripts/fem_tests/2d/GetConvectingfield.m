function [Cx Cy] = GetConvectingfield(xm1,ym1);


Cx=zeros(size(xm1));
Cy=zeros(size(ym1));

U0=1.0;
du=0.1*U0;

alpha=pi/2;
beta=pi/2;

V0=1.0;
dv=-du*alpha/beta;

Cx = U0 + du*sin(alpha*xm1 + beta*ym1);
Cy = V0 + dv*sin(alpha*xm1 + beta*ym1);

% Might need to project lower polynomial order.
% Does it still remain divergence free?

