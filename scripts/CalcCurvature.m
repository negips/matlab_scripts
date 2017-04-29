function [k] = CalcCurvature(x,y);

dydx = gradient(y,x);
d2ydx2 = gradient(dydx,x);

k = abs(d2ydx2)./((1 + dydx.^2).^3/2);

return
