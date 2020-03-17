%     First attempts at modifying the mesh

clear
clc
close all
%
%ifplot = 0;
%casename = 'saab750k';
casename = 'naca0012_2d';

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

OFP = 15.0;
DX  = -5.0;
nsteps = 100;
dx  = DX/nsteps;

xc_old = rea.mesh.xc;

for i=1:nsteps
  xc2 = OFP - rea.mesh.xc; 
  DX2 = dx*blending_erf(xc2);
  OFP = OFP + dx;
  
  rea.mesh.xc = rea.mesh.xc + DX2;
end

CreateVTKMesh(rea.mesh);

xin = 10 - rea.mesh.xc(:);
xin2 = sort(xin);
yin  = blending_erf(xin2);
semilogy(xin2,(yin+0.99e-8)); grid on

xdiff = xc_old - rea.mesh.xc;
[xx ind] = sort(xc_old(:)); 
xdiff = xdiff(ind);

figure
semilogy(xx,(xdiff+1.0e-15)); grid on

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
    DX2 = dx*blending_erf(xgll2);
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


      
