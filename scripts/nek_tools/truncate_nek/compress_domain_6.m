%     First attempts at modifying the mesh

%                 f2
%           x3-----------x2      
%           |            |
%           |            |
%         f3|            |f1
%           |            |
%           |            |
%           x4-----------x1
%                 f4

clear
clc
close all
%
%ifplot = 0;
%casename = 'saab750k';
%casename = 'naca0012_2d';
casename = 'test2d';

rea = Nek_ReadRea(casename);

GLOLD = 1:rea.mesh.nelg;

% Find outflow boundary elements
% Asumming one element has atmost one outflow BC
disp('Finding Outflow Elements')
ELOute = [];      % Element number with Outflow BC
ELOutf = [];      % Outflow BC on face no
nfaces = 2*rea.mesh.ndim;
for e=1:rea.mesh.nelg
  for f=1:nfaces
    cb = rea.mesh.cbc(f,e).bc;
    ifout = strcmpi(cb,'o  ');
    if ifout
      ELOute = [ELOute e];
      ELOutf = [ELOutf f];
      break
    end
  end
end  

nout = length(ELOute);
disp([num2str(nout) ' Outflow Elements.'])

% %
disp('Finding Dirichlet Elements')
ELVe = [];      % Element number with Outflow BC
ELVf = [];      % Outflow BC on face no
nfaces = 2*rea.mesh.ndim;
for e=1:rea.mesh.nelg
  for f=1:nfaces
    cb = rea.mesh.cbc(f,e).bc;
    ifv = strcmpi(cb,'v  ');
    if ifv
      ELVe = [ELVe e];
      ELVf = [ELVf f];
      break
    end
  end
end  

nv = length(ELVe);
disp([num2str(nv) ' Dirichlet Elements.'])

% %
disp('Finding wall/mv Elements')
ELWe = [];      % Element number with Outflow BC
ELWf = [];      % Outflow BC on face no
nfaces = 2*rea.mesh.ndim;
for e=1:rea.mesh.nelg
  for f=1:nfaces
    cb = rea.mesh.cbc(f,e).bc;
    ifW  = strcmpi(cb,'W  ');
    ifmv = strcmpi(cb,'mv ');
    if ifW || ifmv
      ELWe = [ELWe e];
      ELWf = [ELWf f];
      break
    end
  end
end  

nw = length(ELWe);
disp([num2str(nw) ' Wall/mv Elements.'])

% Find nearest Outflow-boundary points.
xoutf = [];
youtf = [];
elout = [];
fcout = [];

for i=1:length(ELOute)
  j=ELOute(i);
  
  elout = [elout; ELOute(i); ELOute(i)];
  fcout = [fcout; ELOutf(i); ELOutf(i)];
 
  if ELOutf(i)==1
    xoutf = [xoutf; rea.mesh.xc(1,j); rea.mesh.xc(2,j)];
    youtf = [youtf; rea.mesh.yc(1,j); rea.mesh.yc(2,j)];
  elseif ELOutf(i)==2
    xoutf = [xoutf; rea.mesh.xc(2,j); rea.mesh.xc(3,j)];
    youtf = [youtf; rea.mesh.yc(2,j); rea.mesh.yc(3,j)];
  elseif ELOutf(i)==3
    xoutf = [xoutf; rea.mesh.xc(3,j); rea.mesh.xc(4,j)];
    youtf = [youtf; rea.mesh.yc(3,j); rea.mesh.yc(4,j)];
  else
    xoutf = [xoutf; rea.mesh.xc(3,j); rea.mesh.xc(4,j)];
    youtf = [youtf; rea.mesh.yc(3,j); rea.mesh.yc(4,j)];
  end        
end

[r c] = size(rea.mesh.xc);
npts = r*c;

dmin_el      = 0*rea.mesh.xc;
dmin_fc      = 0*rea.mesh.xc;
dmin_pt      = 0*rea.mesh.xc;
dmin_ds      = 0*rea.mesh.xc;

for el=1:rea.mesh.nelg
  for i=1:4
    dist = sqrt((xoutf - rea.mesh.xc(i,el)).^2 + (youtf - rea.mesh.yc(i,el)).^2);
    [val ind] = min(dist);
    el2 = elout(ind);
    fc2 = fcout(ind);
    pt2 = mod(ind,2);
    if (pt2==0)
       pt2 = 2;
    end
    dmin_ds(i,el) = val;
    dmin_el(i,el) = el2;
    dmin_fc(i,el) = fc2;
    dmin_pt(i,el) = (fc2-1) + pt2;
    if dmin_pt(i,el)==5
      dmin_pt(i,el)=1;
    end  
%    dmin_pt(i,el) = pt2;
  end
end  
            
OFP = 5.0;
DX  = 1.0;
nsteps = 1000;
dx  = DX/nsteps;

xc_old = rea.mesh.xc;
%
%for i=1:nsteps
%  xc2 = OFP - rea.mesh.xc; 
%  DX2 = (xc2/nsteps).*blending_exp(xc2);
%%  OFP = OFP + dx;
%
%  rea.mesh.xc = rea.mesh.xc + DX2;
%end

xc2 = OFP - rea.mesh.xc;
for el=1:rea.mesh.nelg
  for i=1:4
     el2=dmin_el(i,el);
     pt2=dmin_pt(i,el);
     xc2(i,el) = xc2(pt2,el2);
  end
end  

DX2 = xc2.*blending_exp(dmin_ds);

rea.mesh.xc = rea.mesh.xc + DX2;

CreateVTKMesh(rea.mesh);

xin = dmin_ds(:);
xin2 = sort(xin);
yin  = blending_exp(xin2);

figure
semilogy(xin2,(yin+0.99e-8)); grid on
title('blending function with distance')

xdiff = xc_old - rea.mesh.xc;
[xx ind] = sort(xc_old(:)); 
xdiff = xdiff(ind);

figure
semilogy(xx,abs(xdiff+1.0e-15)); grid on

Nek_WriteRea(rea,0);

%files = {'naca00120.f00002'};
files = {''};
nfiles = length(files);

for ii=1:nfiles

  fname = files{ii};
  if isempty(fname)
    continue
  end  
  [datastruc,data,lr1,nelt,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(fname);
  
  xgll = datastruc(1).data;
  OFP = 15.0;
  DX  = -5.0;
  nsteps = 100;
  dx  = DX/nsteps;
  
  for i=1:nsteps
    xgll2 = OFP - xgll; 
    DX2 = dx*blending_exp(xgll2);
    OFP = OFP + dx;
    
    xgll = xgll + DX2;
  end

  ygll = datastruc(2).data;
  zgll = datastruc(3).data;
  
  Glno = elmap;
  U = datastruc(4).data;
  V = datastruc(5).data;
  W = datastruc(6).data;
  P = datastruc(7).data;
  T = W;
  ifx = 1;
  ifu = 1;
  ifp = 1;
  nps = 0;
  ndim = 3;
  N = 8;
  nel = length(Glno);
  wdsiz = 4;

  fname = sprintf('%s%5.5i','naca00120.f',ii);
  [status] = Nek_WriteFld(ndim,N,nel,xgll,ygll,zgll,U,V,W,P,T,nps,ifx,ifu,ifp,Glno,wdsiz,fname)

end  % ii=1:nfiles 



%====================================================================== 
function blend = blending_erf(x)

  tol=1.0e-8;

  xs = -1.0;
  xe = 0.0;
  ys = erf(xs);
  ye = erf(xe);

  xst = 0;
  xen = 9;
  x = (x-xst)/(xen-xst)*(xe-xs) + xs;

  blend = (ye - erf(x))./(ye - ys);
  ind = find(blend<tol);
  blend(ind) = 0.;

end

%====================================================================== 
function blend = blending_exp(x)

  tol=1.0e-8;

  mu = 0.65;

  blend = exp(-(x/mu).^2);
  ind = find(blend<tol);
  blend(ind) = 0.;

end

      
