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

% % Just plotting
% fig1 = figure(1);
% for i=1:nout
%   el = ELOute(i);
%   Plot2DElement(rea.mesh,el,fig1);
% end  
% for i=1:nv
%   el = ELVe(i);
%   Plot2DElement(rea.mesh,el,fig1);
% end  

% Remove n layers from outflow
ELrm = [];
ELRM = zeros(1,rea.mesh.nelg);

n=11;

Oe = ELOute;
Of = ELOutf;

for il=1:n
  LayerOe{il} = Oe;
  LayerOf{il} = Of;
     
  for e=1:nout
    el = LayerOe{il}(e);
    f  = LayerOf{il}(e);

    ELrm = [ELrm el];         % Remove these elements
    ELRM(el) = 1;

    fopp = mod(f+2,nfaces);
    if fopp==0
      fopp=4;
    end
    el = rea.mesh.cbc(fopp,el).connectsto;      % Element in the next layer
    f  = rea.mesh.cbc(fopp,el).onface;          % Face no of next element

    Oe(e) = el;
    Of(e) = f;
  end
end  


% Remove n layers from outflow
n=3;

Ve = ELVe;
Vf = ELVf;

for il=1:n
  LayerVe{il} = Ve;
  LayerVf{il} = Vf;
     
  for e=1:nv
    el = LayerVe{il}(e);
    f  = LayerVf{il}(e);

    ELrm = [ELrm el];         % Remove these elements
    ELRM(el) = 1;

    fopp = mod(f+2,nfaces);
    if fopp==0
      fopp=4;
    end
    el = rea.mesh.cbc(fopp,el).connectsto;      % Element in the next layer
    f  = rea.mesh.cbc(fopp,el).onface;          % Face no of next element

    Ve(e) = el;
    Vf(e) = f;
  end
end  

ELrmq = unique(ELrm);

nrm = length(ELrmq);

disp([num2str(nrm) ' Elements to be removed'])

csum = cumsum(ELRM);
%plot(csum)

GLNEW = GLOLD - csum;

rmind = find(ELRM>0);
keepind = find(ELRM==0);
GLNEW2 = GLNEW;
GLNEW2(rmind) = [];
GLNEW(rmind) = -1;

% Update Element numbers

rea2 = rea;

for e=1:rea2.mesh.nelg
%  Update Global numbers      
   rea2.mesh.globalno(e) = GLNEW(e);
   
%  Update CBC array   
   for j=1:nfaces
     el = rea2.mesh.cbc(j,e).connectsto;
     if (el>0)
       rea2.mesh.cbc(j,e).connectsto = GLNEW(el);
       if (GLNEW(el)<0)
         rea2.mesh.cbc(j,e).onface = 0;
       end
     end
   end

end   

% Update Curved Element numbers
for e=1:rea2.mesh.Ncurve
  el = rea2.mesh.curveieg(e);
  rea2.mesh.curveieg(e) = GLNEW(el);
end

% Update Outflow Boundary conditions
for i=1:nout
  el = Oe(i);
  f  = Of(i);

  rea2.mesh.cbc(f,el).bc = 'O  ';
  rea2.mesh.cbc(f,el).connectsto = 0;
  rea2.mesh.cbc(f,el).onface = 0;
end  

% Update Outflow Boundary conditions
for i=1:nv
  el = Ve(i);
  f  = Vf(i);

  rea2.mesh.cbc(f,el).bc = 'v  ';
  rea2.mesh.cbc(f,el).connectsto = 0;
  rea2.mesh.cbc(f,el).onface = 0;
end  

rea2.mesh.globalno = rea2.mesh.globalno(keepind);
rea2.mesh.groupno  = rea2.mesh.groupno(keepind);
rea2.mesh.xc = rea2.mesh.xc(:,keepind);
rea2.mesh.yc = rea2.mesh.yc(:,keepind);

rea2.mesh.cbc = rea2.mesh.cbc(:,keepind);

ind = find(rea2.mesh.curveieg>0);
rea2.mesh.Ncurve=length(ind);
rea2.mesh.curveieg=rea2.mesh.curveieg(ind);
rea2.mesh.curveedge=rea2.mesh.curveedge(ind);
rea2.mesh.curveparams=rea2.mesh.curveparams(:,ind);
rea2.mesh.curvetype=rea2.mesh.curvetype(ind);

nelg = length(keepind);
rea2.mesh.nelg = nelg;

CheckConnectivity2D(rea2.mesh)


% 
CreateVTKMesh(rea2.mesh);
Nek_WriteRea(rea2,0);
% 
fid = fopen('truncated.numbering', 'w');

for i=1:rea.mesh.nelg
  fprintf(fid,'%i %i\n', GLOLD(i), GLNEW(i));
end
fclose(fid)
% 
clearvars -except GLOLD GLNEW

files = {'bsenaca00120.f00001'};
%files = {''};
nfiles = length(files);

for ii=1:nfiles

  fname = files{ii};
  if isempty(fname)
    continue
  end  
  [datastruc,data,lr1,nelt,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(fname);
  
  Elmap2 = 0*GLOLD;
  for i=1:nelt
    no = elmap(i);
    Elmap2(i) = GLNEW(i);
  end  
  
  IndT = find(Elmap2>0);
  GLNEW2 = Elmap2(IndT);
  
  % Remove unwanted elements
  nf = length(datastruc);
  for i=1:nf
    datastruc(i).data = datastruc(i).data(:,IndT);
  end  
  
  [GLNEW3 I] = sort(GLNEW2);
  
  % reorder data
  for i=1:nf
    datastruc(i).data = datastruc(i).data(:,I);
  end 
  
  Glno = GLNEW3;
  xgll = datastruc(1).data;
  ygll = datastruc(2).data;
  zgll = 0*datastruc(2).data;
  U = datastruc(3).data;
  V = datastruc(4).data;
  W = 0*datastruc(2).data;
  P = datastruc(5).data;
  T = W;
  ifx = 1;
  ifu = 1;
  ifp = 1;
  nps = 0;
  ndim = 2;
  N = 10;
  nel = length(Glno);
  wdsiz = 8;

  fname = sprintf('%s%5.5i','naca00120.f',ii);
  [status] = Nek_WriteFld(ndim,N,nel,xgll,ygll,zgll,U,V,W,P,T,nps,ifx,ifu,ifp,Glno,wdsiz,fname)

end  % ii=1:nfiles 



