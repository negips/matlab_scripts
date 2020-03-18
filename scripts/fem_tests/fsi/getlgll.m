function [xl yl] = getlgll(Nx,Ny,xc,yc);
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

