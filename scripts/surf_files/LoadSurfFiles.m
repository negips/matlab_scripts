function [fout tout] = LoadSurfFiles(fol)

command = ['ls ' fol '/surf_data01_*'];
[status,cmdout] = system(command);

ll = length(cmdout);

l1 = length(fol);
samplefilename = 'surf_data01_0.10000000E+01.bin';
l2 = length(samplefilename);
prefix = 'surf_data01_';
l3 = length(prefix);

flen = l1+l2+1;         % add 1 for '/'

inds = strfind(cmdout,fol);

nfiles = length(inds);

tout = [];
for i = 1:nfiles
  ind1=inds(i);
  ind2=ind1+flen-1;
  file = cmdout(ind1:ind2); 
  allfiles{i}=file;
  time = file(l1+l3+2:flen-4);
  tout(i)=str2num(time);
end

[tout ind] = sort(tout);

fout = allfiles(ind);

return

