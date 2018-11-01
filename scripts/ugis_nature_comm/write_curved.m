%% Creating curved side data (mid point)

clear
clc
close all

surf = importdata('surf_data.dat');

data = surf.data;
N=7;
lx1=8;

npts = length(data(:,1));
nels = npts/lx1

x=zeros(lx1,nels);
y=zeros(lx1,nels);
face=zeros(nels);
ieg=zeros(nels,1);
ifcirc=zeros(nels,1);
ifround=zeros(nels,1);
ncirc=0;
nround=0;   % for rounded edges

xcut = 0.4998999;
xcut_round = 1.5;

for iel=1:nels
  for i=1:lx1
    ix = (iel-1)*lx1 + i;
    x(i,iel)=data(ix,1);
    y(i,iel)=data(ix,2);
    glno = data(ix,7);  
    f=data(ix,8);
  end
  face(iel)=f;
  ieg(iel)=glno;
  if max(x(:,iel)>xcut)
    ifcirc(iel) = 0;
  else
    ncirc = ncirc + 1;
    ifcirc(iel) = 1;
  end

% Rounded end
  if max(x(:,iel)>xcut_round)
    nround = nround + 1;
    ifround(iel) = 1;
  else
    ifround(iel) = 0;
  end

end  

ind = 1:nels;
circind=ind(find(ifcirc));

xcirc = x(:,circind);
ycirc = y(:,circind);

endx=xcirc([1 lx1],:);
endy=ycirc([1 lx1],:);

%endtheta = atan2(endx/a,endy/b);

xc = mean(endx,1);
yc = mean(endy,1);

pgll = plot(x,y,'.-b'); hold on
pend = plot(endx(:),endy(:), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plmid= plot(xc(:),yc(:), 'gs', 'MarkerSize', 12);

a=0.5;
b=0.5;

theta = atan2(yc/b,xc/a);
xmid = a*cos(theta);
ymid = b*sin(theta);

pcmid = plot(xmid,ymid, 'md', 'MarkerSize', 12);

fcirc=face(circind);
iegcirc=ieg(circind);


% Rounded edge
roundind=ind(find(ifround));
xround = x(:,roundind);
yround = y(:,roundind);

endx2=xround([1 lx1],:);
endy2=yround([1 lx1],:);

xc2 = mean(endx2,1);
yc2 = mean(endy2,1);

pend2 = plot(endx2(:),endy2(:), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plmid2= plot(xc2(:),yc2(:), 'gs', 'MarkerSize', 12);

a=0.01;
b=0.01;

X0=1.499899;
Y0=0.;

theta = atan2((yc2-Y0)/b,(xc2-X0)/a);
xmid2 = X0 + a*cos(theta);
ymid2 = Y0 + b*sin(theta);

pcmid2 = plot(xmid2,ymid2, 'md', 'MarkerSize', 12);

fround=face(roundind);
ieground=ieg(roundind);

legend([pgll(1), pend, plmid, pcmid], {'GLL', 'End-Pts', 'Avg.', 'Mid-Pts'})

% Modify for rounded end

% cc = "{0:6d} Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE \n"
% c1 = "{0:2d}{1:6d}"
%if tot_num_cells < 1000:
%    c1 = "{0:3d}{1:3d}"
%else:
%    c1 = "{0:2d}{1:6d}"
%end
% c2 = "{0:14.6e}{1:14.6e}{2:14.6e}{3:14.6e}{4:14.6e} {5:s}\n"

if nels<1000
   fmt='%2d%6d%14.6e%14.6e%14.6e%14.6e%14.6e %s\n';
else
   fmt='%3d%3d%14.6e%14.6e%14.6e%14.6e%14.6e %s\n';
end

fid = fopen('curved.out', 'w');
fprintf(fid,'%6d Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n',(ncirc+nround));
for i=1:ncirc
  fprintf(fid,fmt, fcirc(i), iegcirc(i),xmid(i),ymid(i), 0.0, 0.0, 0.0, 'm');
end  
for i=1:nround
  fprintf(fid,fmt, fround(i), ieground(i),xmid2(i),ymid2(i), 0.0, 0.0, 0.0, 'm');
end  

fclose(fid);








