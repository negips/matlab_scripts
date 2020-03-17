close all
clear all
clc
double precision;

%Initialize nponits
npoints=0;

%Generate interpolating mesh
%Region 1
% xx=-1:0.01:1;
% yy=-1:0.01:1;

xx=-1:0.005:1;
yy=-1:0.005:1;

for i=1:length(xx)
    for j=1:length(yy)
        npoints=npoints+1;
        x_pts(npoints,1)=xx(i);
        y_pts(npoints,1)=yy(j);
    end
end

npoints

%Save x data points in Fortran binary format
fid=fopen('./ZSTAT/x.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=x_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%Save y data points in Fortran binary format
fid=fopen('./ZSTAT/y.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=y_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

figure(1)
plot(x_pts,y_pts,'.b'); hold on;

axis equal

save('int_mesh.mat');