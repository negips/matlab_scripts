function xl = getlgll1D(Nx,xs,xe)
%    Element mapping in 1D

x = lglnodes(Nx);
x = x(end:-1:1);

jx = zeros(length(x),2);
dl = x(end) - x(1);
for i=1:length(x);
    jx(i,2) = (x(i) -x(1))/dl;
    jx(i,1) = (x(end) - x(i))/dl; 
end

%    mxm(dx,lx1,xcb,2,w,2)

%xcb = [[-1; 2] [0; 4]];
xcb = zeros(2,1);
xcb(1) = xs;
xcb(2) = xe;
xl = jx*xcb;

