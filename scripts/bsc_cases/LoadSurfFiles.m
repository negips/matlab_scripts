function [fout tout] = LoadSurfFiles(fol)

command = ['ls ' fol '/surf_data01_*'];
[status,cmdout] = system(command);

ll = length(cmdout);

l1 = length(fol);
samplefilename = 'surf_data01_0.10000000E+01.bin';
l2 = length(samplefilename);

flen = l1+l2;

inds = strfind(cmdout,fol);

nfiles = length(inds);
prefix = 'surf_data01_';
l3 = length(prefix)

tout = [];
for i = 1:nfiles
  ind1=inds(i)+l1+1;
  ind2=ind1+l2;
  file = cmdout(ind1:ind2); 
  allfiles{i}=file;
  tout(i)=str2num(file(l3+1:end-5));
  file
      
end

[tout ind] = sort(tout);

fout = allfiles(ind);

return

