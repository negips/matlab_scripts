function [xl yl] = getlgll_ho(Nx,Ny,xc,yc);
%    testing element mapping

addpath '../../'

ind1 = [1 2];
ind2 = [4 3];

x = lglnodes(Nx);
x = x(end:-1:1);

jx = zeros(length(x),2);
dl = x(end) - x(1);
for i=1:length(x);
    jx(i,2) = (x(i) -x(1))/dl;
    jx(i,1) = (x(end) - x(i))/dl; 
end

jxt = transpose(jx);

y = lglnodes(Ny);
y = y(end:-1:1);

jy = zeros(length(y),2);
dl = y(end) - y(1);
for i=1:length(y);
    jy(i,2) = (y(i) -y(1))/dl;
    jy(i,1) = (y(end) - y(i))/dl; 
end

jyt = transpose(jy);

%    mxm(dx,lx1,xcb,2,w,2)

%xcb = [[-1; 2] [0; 4]];
xcb = zeros(2,2);
xcb(:,1) = xc(ind1);
xcb(:,2) = xc(ind2);
w = jx*xcb;
xl = w*jyt;

ycb = zeros(2,2);
ycb(:,1) = yc(ind1);
ycb(:,2) = yc(ind2);
w = jx*ycb;
yl = w*jyt;


ifexp   = 1;
ifparab = 0;
ifcubic = 0;

if (ifexp)
  % Rotate points with radial exponential decay
  
  % Axis of rotation
  x0=0.5;
  y0=0.5;
  
  xl2 = xl(:)-x0;
  yl2 = yl(:)-y0;
  
  theta = 50*pi/180;
  Rot   = [cos(theta) -sin(theta);...
           sin(theta)  cos(theta)];
  
  coords = Rot*[xl2'; yl2'];
  xl3 = reshape(coords(1,:),Nx+1,Ny+1) + x0;
  yl3 = reshape(coords(2,:),Nx+1,Ny+1) + y0;
  
  dx  = xl3 - xl;
  dy  = yl3 - yl;
  
  dist  = sqrt((xl-x0).^2 + (yl-y0).^2);
  
  mu    = 0.2;
  decay = exp(-(dist/mu).^4);
  dx    = dx.*decay;
  dy    = dy.*decay;
  
  xl    = xl + dx;
  yl    = yl + dy;

elseif (ifparab)

  % 2nd order displacement for x/y points.
  
  % Maximum displacement 
  dx0=0.2;
  dy0=0.2;
  
  dx = dx0*((yl-0.5)/0.5).^2;
  dy = dy0*((xl-0.5)/0.5).^2;
  
  xl  = xl + dx;
  yl  = yl + dy;

elseif (ifcubic)

  % 3rd order displacement for x/y points.
  
  % Maximum displacement 
  dx0=0.1;
  dy0=0.1;
  
  dx = dx0*(abs(yl-0.5)/0.5).^3;
  dy = dy0*(abs(xl-0.5)/0.5).^3;
  
  xl  = xl + dx;
  yl  = yl + dy;

end














