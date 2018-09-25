% Modify template .rpl file to add spline to points

clear
clc

fid=fopen('add_curve.rpl');
fid2 = fopen('add_new_curve2.rpl','w');

tline = fgetl(fid);
npts = 1000; 
while ischar(tline)
  ind = findstr(tline, '{pnt_fine0'); 
  if isempty(ind)
    fprintf(fid2,'%s\n',tline);
  else
    tline2 = tline(1:ind);  
    for i=0:npts-1
       tline2 = [tline2 'pnt_v2' num2str(i) ' '];
    end
    tline2 = [tline2 '}'];
    fprintf(fid2,'%s\n', tline2);
  end
  tline = fgetl(fid);
end

fclose all; 
          

