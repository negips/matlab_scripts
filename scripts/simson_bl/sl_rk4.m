function [xsl ysl] = sl_rk4(x,y,u,v,h,x0,y0)

%  h;                                             % step size
  
%  xk_1 = interp2(x,y,u,x0,y0);
%  yk_1 = interp2(x,y,v,x0,y0);

  x1 = x0+0.5*h;
  y1 = y0+0.5*h; 

%  xk_2 = interp2(x,y,u,x1,y1);
%  yk_2 = interp2(x,y,v,x1,y1);

  x2 = x0+0.5*h;
  y2 = y0+0.5*h;

%  xk_3 = interp2(x,y,u,x2,y2);
%  yk_3 = interp2(x,y,v,x2,y2);
     
  x3 = x0+h;
  y3 = y0+h;

%  xk_4 = interp2(x,y,u,x3,y3);
%  yk_4 = interp2(x,y,v,x3,y3);

  xkn = interp2(x,y,u,[x0 x1 x2 x3], [y0 y1 y2 y3]);
  ykn = interp2(x,y,v,[x0 x1 x2 x3], [y0 y1 y2 y3]);
 
  xsl = x0 + (1/6)*(xkn(1) + 2*xkn(2) + 2*xkn(3) + xkn(4))*h;
  ysl = y0 + (1/6)*(ykn(1) + 2*ykn(2) + 2*ykn(3) + ykn(4))*h;

%  xsl = x0 + (1/6)*(xk_1+2*xk_2+2*xk_3+xk_4)*h;  % main equation
%  ysl = y0 + (1/6)*(yk_1+2*yk_2+2*yk_3+yk_4)*h;  % main equation


return
