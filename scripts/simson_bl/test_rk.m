% testing rk_implementation

x=0:1000;
y=0:100;

[X Y] = meshgrid(x,y);
U = ones(size(X));
V = 0.015*U;

x_seed = 1; 
nseed = 1;
dr = 1e-2;
maxverts=500;
y_seed = 0.5;

x_sl = [];
y_sl = [];
y_sl2 = [];
xsl_max = 500;
ysl_max = 50;

for iseed = 1:nseed
  xsl=x_seed;
  ysl=y_seed(iseed);
  xsl_arr = [xsl];
  ysl_arr = [ysl];
  npts=0;    
  while ((xsl<xsl_max) && (ysl<ysl_max))
    clc  
    npts=npts+1
    iseed  
    [xsl ysl] = sl_rk4(X,Y,U,V,dr,xsl,ysl);
    xsl_arr = [xsl_arr; xsl];
    ysl_arr = [ysl_arr; ysl];
    if npts>maxverts
      break
    end  
  end
  slines(iseed).x=xsl_arr;
  slines(iseed).y=ysl_arr;
  slines(iseed).xw = x;
end     

figure(20)
plot(xsl_arr,ysl_arr)

