% Create 3D Channel box file

clear
clc

format long

load 'boxpar.mat'

ndim=3;
nfld=1;
fname='moser395';
reaname='turbChannel'
%nelx=32;
%nely=16;
%nelz=16;
bcs='PPWWPP';


fmt1='%12.9f ';
fmt2='%12.9f \n';
fmt3='%6.3f ';
fmt4='%6.3f \n';

pi=4*atan(1);

boxfile=[fname '.box2'];
reafile=[reaname '.rea'];

fid=fopen(boxfile, 'w');

fprintf(fid,'%s\n',reafile);
fprintf(fid,'%d\t\t\t Spatial Dimension\n',ndim);
fprintf(fid,'%d\t\t\t Number of fields\n',nfld);
fprintf(fid,'#\t\t\t Comments\n');
fprintf(fid,'#\n');
fprintf(fid,'#\n');
fprintf(fid,'#===========================================\n');
fprintf(fid,'#\n');
fprintf(fid,'Box 1\t\t\t\t\t\t Channel\n');
fprintf(fid,'%d %d %d\t\t\t\t\t Nelx Nely Nelz\n',-nelx,nely,-nelz);


%% Xpoints
x0=0.0;
x1=Lx;
xr=1.0;
fprintf(fid,'%4f %4f %4f \t\t\t x0 x1 ratio\n',x0,x1,xr);

%% Y points
cnt1=0;
y0=-1;
y1=Ly-1;
yrange=y1-y0;
ymid=(y0+y1)/2;

for i=0:nely
    
    if cnt1>90
        fmt=fmt4;
        cnt1=0;
    else
        fmt=fmt3;
    end
    
%    y_pts(i+1)=ymid+yrange/2*cos(pi/(nely)*i-pi);
    fprintf(fid,fmt,ypts(i+1));
    
    cnt1=cnt1+7;
    
end

if strcmp(fmt,fmt3)
    fprintf(fid,'\n');
end

disp(['Y points:'])
ypts'

%% Z points

z0=0.0;
z1=Lz;
zr=1.0;
fprintf(fid,'%5f %5f %5f \t\t\t z0 z1 ratio\n',z0,z1,zr);
%% Boundary Conditions

for i=1:2*ndim-1
    fprintf(fid,'%s  ,',bcs(i));   
end

fprintf(fid,'%s  ',bcs(i+1));   
fprintf(fid,'\t\t\t\t BC''s: (x,y,z)');

fclose(fid);





