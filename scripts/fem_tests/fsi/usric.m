function uic = usric(x,y)

     pi = 3.14159265;
     n=5;
     m=5;
%     uic = 1/4*(4*y-y^2)*sin(2*pi*x);
%     uic = sin(2*pi*x/3);
     uic = 0.0 + 0*sin(n*2*pi*x)*sin(m*2*pi*y);

     % Build a Gaussian initial condition
%     u0 = 1.0;
%     sigmax = 0.1;
%     sigmay = 0.1;
%     x0 = 0.5;
%     y0 = 0.5;
%     uic = u0*exp( -(((x-x0)/sigmax).^2) - (( (y-y0)/sigmay).^2));

return
