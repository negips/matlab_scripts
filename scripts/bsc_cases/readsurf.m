function [sdata sintegrals tstamps sno lx1 selt maxtsaves x y timeout hdr] = readsurf(fname,ifhdr)

sdata=[];
sintegrals = [];
tstamps = [];
sno = [];
lx1 = [];
selt = [];
maxtsaves = [];
x = [];
y = [];
timeout = [];
hdr = [];

en='le';
[fid,message] = fopen(fname,'r',['ieee-' en]);

if ifhdr

  nread=132;
  hdr=fread(fid,132,'*char');
  hdr = transpose(hdr);
  dummy = sscanf(hdr,['%s%d%d%d%d%s%d%f'] );
  hash = char(dummy(1));
  precision = dummy(2);
  sno = dummy(3);
  lx1 = dummy(4);
  selt = dummy(5);
  ifx = char(dummy(6)); ifx = 1;    
  ifs = char(dummy(7)); ifs = 1;
  ns1 = 10*str2num(char(dummy(8))); 
  ns2 = 1*str2num(char(dummy(9)));
  ns = ns1+ns2;
  maxtsaves = dummy(10);
  timeout = dummy(11);

  if precision == 8 
     fp='*float64';
  elseif precision == 4
     fp='*float32';
  else
     disp(['Unknow precision:' num2str(precision)])
     fclose(fid);
     return  
  end              


else

  %% New temporary header
  hdr='No header';
  fp='*float64';
  
  nread=1;
  
  precision = fread(fid,nread,fp);
  precision = precision(1);
  
  % surface no
  sno = fread(fid,nread,fp);
  sno = sno(1);
  
  % Polynomial order
  lx1 = fread(fid,nread,fp);
  lx1 = lx1(1);
  
  % No of surface elements
  selt = fread(fid,nread,fp);
  selt = selt(1);
  
  % No of variables
  ns = fread(fid,nread,fp);
  ns = ns(1);
  
  % No of time history points.
  maxtsaves = fread(fid,nread,fp);
  maxtsaves = maxtsaves(1);
  
  % Output time.
  timeout = fread(fid,nread,fp);
  timeout = timeout(1);
  
  ifx=1;
  ifs=1;

end

displayhdr = 1;
if displayhdr
  disp(['Precision: ' num2str(precision)])    
  disp(['Surf No: ' num2str(sno)])    
  disp(['Selt: ' num2str(selt)])    
  disp(['lx1: ' num2str(lx1)])    
  disp(['Num. of variables: ' num2str(ns)])    
  disp(['Max Time saves: ' num2str(maxtsaves)])    
  disp(['Output time: ' num2str(timeout)])
  disp('---------------------')
end 

%dbstop in readsurf at 88

%% Read Time stamps 
tstamps=zeros(maxtsaves,1);
for i=1:maxtsaves
  dtmp=fread(fid,1,fp);
  tstamps(i) = dtmp(1);
end    

%---------------------------------------- 
ndim=3;
sintegrals = zeros(maxtsaves,ndim,ns);
for i=1:ns
  for j=1:ndim
    for k=1:maxtsaves  
      tmp=fread(fid,1,fp);
      sintegrals(k,j,i)=tmp(1);
    end  
  end
end
%---------------------------------------- 

%dbstop in readsurf at 55

x=zeros(lx1,selt);
y=zeros(lx1,selt);

if (ifx)
  for i=1:selt
    for j=1:lx1  
      tmp=fread(fid,1,fp);
      x(j,i)=tmp(1);
    end
    
    for j=1:lx1
       tmp=fread(fid,1,fp);
       y(j,i)=tmp(1);
    end      
  end
end

%---------------------------------------- 

sdata = [];
scdata=zeros(lx1,selt,maxtsaves);
if (ifs)
  for sc=1:ns
    for it=1:maxtsaves
      for i=1:selt
        for j=1:lx1
          tmp=fread(fid,1,fp);
          scdata(j,i,it)=tmp(1);
        end    
      end
    end
    sdata(sc).data=scdata;
  end
end


fclose(fid);

% testing
% h1=figure;
% hold on
% for i=1:selt
% %     scatter(x(:,i),sdata(1).data(:,i,1), '.')
% %     scatter(x(:,i),sdata(1).data(:,i,maxtsaves), 'd')
%      scatter(sdata(1).data(:,i,maxtsaves),sdata(3).data(:,i,maxtsaves), 'd')
% 
% end
% set(gca,'Ydir','reverse')

return

