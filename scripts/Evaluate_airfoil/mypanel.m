function [u,cp]=mypanel(x,y,alpha)
%  [u,cp]=mypanel(x,y,alpha)
%
%  Panel Code in MATLAB
%
%  Based onCoded by L. sankar, April 1997
%
%  Input:
%  x = x-xoordinates of geometry
%  y = y-coordinate of geometry
%  alpha = angle of attack in degree
%

n=size(x,1)-1;
a=zeros(n+1,n+1);
%
%compute arclength
%
dx=diff(x);
dy=diff(y);
ds=sqrt(dx.*dx+dy.*dy);
s=[0;cumsum(ds)];

sint=dy./ds;
cost=dx./ds;

xmid=(x(1:end-1)+x(2:end))*0.5;
ymid=(y(1:end-1)+y(2:end))*0.5;

%
% Assemble the Influence Coefficient Matrix A
%
for j = 1:n
    a(j,n+1) = 1.0;
    for i = 1:n
        if i == j
            a(i,i) = ds(i)/(2.*pi) *(log(0.5*ds(i)) - 1.0);
        else
            t1  = x(i) - xmid(j);
            t2  = y(i) - ymid(j);
            t3  = x(i+1) - xmid(j);
            t7  = y(i+1) - ymid(j);
            t4  = t1 * cost(i) + t2 * sint(i);
            t5  = t3 * cost(i) + t7 * sint(i);
            t6  = t2 * cost(i) - t1 * sint(i);
            t11  = t5 * log(t5*t5+t6*t6) - t4 * log(t4*t4+t6*t6);
            t12  = atan2(t6,t4)-atan2(t6,t5);
            a(j,i) = (0.5 * t11-t5+t4+t6*t12)/(2.*pi);
        end
    end
    a(n+1,1) = 1.0;
    a(n+1,n) = 1.0;
end
%
% Assemble the Right hand Side of the Matrix system
%
rhs=ymid * cos(alpha* pi /180) - xmid* sin(alpha* pi /180);
rhs=[rhs;0];
%
% Solve the syetm of equations
% In MATLAB this is easy!
%
gamma = a\rhs;
u=abs(gamma(1:end-1));
ds1=(ds(1:end-1)+ds(2:end))*0.5;
s1=[ds(1)/2;cumsum(ds1)+ds(1)/2];
u=interp1(s1,u,s,'pchip');
u=myfilter(s,u,20);
cp=1-u.*u;


end
