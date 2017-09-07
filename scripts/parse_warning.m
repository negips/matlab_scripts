%     Parse the warning file to get coordinates.

clear
clc

fid = fopen('warning.out','r');

npts=0;
while ~feof(fid)
  tline=fgetl(fid);
  A = strread(tline,'%s', 'delimiter', ' ');
  if length(A)==10
    npts=npts+1;
    x(npts)=str2num(A{8});    
    y(npts)=str2num(A{9});    
    z(npts)=str2num(A{10});    
  end
end
fclose(fid)
    
