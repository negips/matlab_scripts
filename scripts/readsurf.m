function [sdata tstamps sno lx1 selt x y timeout hdr] = readsurf(fname)

sdata=[];
sno = -1;

en='le';
[fid,message] = fopen(fname,'r',['ieee-' en]);
reclen = fread(fid,1,'int32');
hdr = fread(fid,reclen,'*char')';
reclen=fread(fid,1,'int32');

% 1002 format (A1,2X,I1,A8,2X,I3.3,2X,I2.2,2X,I5.5,2X,A2,I2.2,2X,I3.3,2X,E14.6E3)

%% parse header
ind=1;
rl=1;
dummy = hdr(ind:ind+rl);      % #
ind=ind+rl+2;
rl=1;
precision = str2double(hdr(ind:ind+rl));
ind=ind+rl+2;
rl=8;
dummy=hdr(ind:ind+rl);        % Surf no
ind=ind+rl+2;
rl=3;
sno = str2double(hdr(ind:ind+rl));
ind=ind+rl+2;
rl=2;
lx1 = str2double(hdr(ind:ind+rl));
ind=ind+rl+2;
rl=5;
selt = str2double(hdr(ind:ind+rl));
ind=ind+rl+2;
rl=4;
flds=hdr(ind:ind+rl);
ifx=0;
ifs=0;
if strfind(flds,'X')
     ifx=1;
end
if strfind(flds,'S')
     ifs=1;
     ns=str2double(flds(3:4));
end
ind=ind+rl+2;
rl=3;
maxtsaves=str2double(hdr(ind:ind+rl));
ind=ind+rl+2;
rl=14;
timeout=str2double(hdr(ind:ind+rl));

if precision==8
     fp='*float64';
else precision==4
     fp='*float32';
end
     

%% Read Time stamps 
tstamps=zeros(maxtsaves,1);
reclen=fread(fid,1,'int32');
tstamps=fread(fid,maxtsaves,fp);
reclen=fread(fid,1,'int32');


x=zeros(lx1,selt);
y=zeros(lx1,selt);


if (ifx)
     for i=1:selt
          reclen=fread(fid,1,'int32');
          x(:,i)=fread(fid,lx1,fp);
          y(:,i)=fread(fid,lx1,fp);
          reclen=fread(fid,1,'int32');
     end
end


sdata = [];
scdata=zeros(lx1,selt,maxtsaves);

if (ifs)
     for sc=1:ns
          for it=1:maxtsaves
               for i=1:selt
                    reclen=fread(fid,1,'int32');
                    tmp=fread(fid,lx1,fp);
                    scdata(:,i,it)=tmp;
                    reclen=fread(fid,1,'int32');
               end
          end
          sdata(sc).data=scdata;
     end
end


fclose(fid);

% testing
h1=figure;
hold on
for i=1:selt
     scatter(x(:,i),sdata(1).data(:,i,1), '.')
     scatter(x(:,i),sdata(1).data(:,i,50), 'd')

end

end

