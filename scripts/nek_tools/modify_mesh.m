%     First attempts at modifying the mesh

clear
clc
close all

casename = 'saab_wing2d';

rea = Nek_ReadRea(casename);

n=rea.mesh.nelg;
ndim=rea.mesh.ndim;
cmap = jet(n);

%figure(1);
%for i=1:n
%  fill(rea.mesh.xc(:,i),rea.mesh.yc(:,i),cmap(i,:)); hold on
%end
%colorbar

% Calculate element mid points.
% Just for plotting for now
xmid=zeros(n,1);
ymid=zeros(n,1);
GLNO=zeros(n,1);
for i=1:n
  xmid(i) = mean(rea.mesh.xc(:,i));
  ymid(i) = mean(rea.mesh.yc(:,i));
  GLNO(i) = rea.mesh.globalno(i);
end  


% Lets do a 2D refine first.

% Find aspect ratios of elements

ARmax = zeros(n,1);
ref_el= zeros(n,1);
ar = zeros(4,n);
rs = zeros(4,n);

arcut=10;    % cut off for refinement/coarsening
for i=1:n
  xt = [rea.mesh.xc(:,i); rea.mesh.xc(1,i)];
  yt = [rea.mesh.yc(:,i); rea.mesh.yc(1,i)];
  dx=diff(xt);
  dy=diff(yt);
  r = sqrt(dx.^2 + dy.^2);
  rs(:,i) = r;
  ar(1,i) = max(r(1)/r(2),r(2)/r(1));
  ar(2,i) = max(r(2)/r(3),r(3)/r(2));
  ar(3,i) = max(r(3)/r(4),r(4)/r(3));
  ar(4,i) = max(r(4)/r(1),r(1)/r(4));
  ARmax(i) = max(ar(:,i));
  xt = rea.mesh.xc(:,i);
  yt = rea.mesh.yc(:,i);
  rad = min(sqrt(xt.^2 + yt.^2));
%  if (rad>3 && ARmax(i)>arcut)
  if (ARmax(i)>arcut)    
    ref_el(i)=1;
  end  
 
end


% find no of neighbors that can be refined 
nfaces=2*ndim;
nbors = zeros(n,1);
for i=1:n
  nrefine=0;    
  for j=1:nfaces
    bc = rea.mesh.cbc(j,i).bc;
    if (strcmpi(bc,'E  '))    
      ieg = rea.mesh.cbc(j,i).connectsto;
      if (ref_el(ieg))
        nrefine=nrefine+1;
      end
    end  
  end
  if nrefine==0
%   If the neighbors can't be modified then the element can't be modified either          
    ref_el(i)=0;
  else
    nbors(i,1)=nrefine;
  end
end

% Find wall/Moving wall layer
wall_el = [];
wall_f  = [];
offwall_el = [];        % Element no of adjoining non-wall element.
ofwall_f   = [];        % face joining this element with wall element
nfaces=2*ndim;
nwall = 0;
for i=1:n
  for j=1:nfaces
    ifwall = 0;    
    cb = rea.mesh.cbc(j,i).bc;
    if (strcmpi(cb,'w  ') || strcmpi(cb,'mv '))
      ifwall = 1;
      wallface = j;
      nwall = nwall+1;
      wall_el(nwall) = i;
      wall_f(nwall) = j;

%     Find adjoining elements without wall-bc
%     These effectively form the layer from the airfoil to the outflow
%      if j==1 || j==2
%        j2=[3 4];
%      else
%        j2=[1 2];
%      end  

      j2=[j+1 j-1];
      if j2(1)>2^ndim
        j2(1)=1;
      end
      if (j2(2)==0)
        j2(2)=4;
      end  
%      offwall_el(nwall) = rea.mesh.cbc(j2,i).connectsto;

      for jj=j2      
        elc = rea.mesh.cbc(jj,i).connectsto;
        haswall = 0;
        for j3=1:nfaces
          cb2 = rea.mesh.cbc(j3,elc).bc;
          if (strcmpi(cb2,'w  ') || strcmpi(cb2,'mv '))
             haswall=1;
          end
        end
        if (~haswall)  
          offwall_el = [offwall_el elc];
        end  
      end   %

    end

  end
end  

cmap = hot(nwall);
figure(2)
for i=1:nwall
  gno = wall_el(i);  
  xt=rea.mesh.xc(:,gno);
  yt=rea.mesh.yc(:,gno);
  fill(xt,yt,cmap(i,:)); hold on
end

offwall_el = unique(offwall_el);
l1 = length(offwall_el);
cmap2 = parula(l1);
for i=1:l1
  gno = offwall_el(i);  
  xt=rea.mesh.xc(:,gno);
  yt=rea.mesh.yc(:,gno);
  fill(xt,yt,cmap2(i,:)); hold on
end








