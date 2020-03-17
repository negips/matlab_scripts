%----------------------------------------------------------------------
%write initial conditions from RANS solution
%----------------------------------------------------------------------

function [status] = Nek_WriteFld(ndim,N,nel,xgll,ygll,zgll,U,V,W,P,T,nps,ifx,ifu,ifp,Glno,wdsz,fname)
      
  lx1 = N+1;

% write header
%  wdsz=8;
  lr1=ones(3,1);
  for i=1:ndim
    lr1(i)=lx1;
  end
  time=0.;
  istep=0.;
  fid=1;    % file no
  nf=1;     % no of files
% Modify fields part of the header
  flds = [];
  if ifx
    flds = [flds 'X'];
  end
  if ifu
    flds = [flds 'U'];
  end
  if ifp
    flds = [flds 'P'];
  end
  if nps>=1
    flds = [flds 'T'];
  end
  if nps>1
    flds = [flds 'S' num2str(nps-1)];
  end  

  hdr = sprintf('#std %1i %2i %2i %2i %10i %10i %20.13E %9i %6i %6i %s\n',...
                 wdsz,lr1(1),lr1(2),lr1(3),nel,nel,time,istep,fid,nf,flds);
  hdr(end+1:132) = ' ';

% word size (double/single precision)
  if (wdsz == 4)
      realtype = 'float32';
  else
      realtype = 'float64';
  end
 
% Open file
  endian= 'le';
  [fid,message] = fopen(fname,'w+',['ieee-' endian]);
  
  if fid == -1
    disp(message) 
    return 
  end
 
  disp(['Writing file: ', fname])

% Write Header  
  fwrite(fid,hdr,'char');

% Endian tag
  etag = 6.54321;
  fwrite(fid,etag,'float32');
  
% Global Element map
  fwrite(fid,Glno,'int32');

% Coordinates
  if ifx
    for ie=Glno
      fwrite(fid,xgll(:,ie),realtype);
      fwrite(fid,ygll(:,ie),realtype);
      if (ndim==3)
        fwrite(fid,zgll(:,ie),realtype);
      end  
    end
  end  

% Velocities  
  if ifu
    for ie=Glno
      fwrite(fid,U(:,ie),realtype);
      fwrite(fid,V(:,ie),realtype);
      if (ndim==3)
        fwrite(fid,W(:,ie),realtype);
      end  
    end
  end  

% Pressure  
  if ifp
    for ie=Glno
      fwrite(fid,P(:,ie),realtype);
    end
  end  

% Temperature  
  for ips=1:nps
    for ie=Glno
      fwrite(fid,T(:,ie,ips),realtype);
    end
  end  

  status = fclose(fid);

end         % function

%----------------------------------------------------------------------
% % element size
% lr1 = [str2double(header(8:9))
%        str2double(header(11:12))
%        str2double(header(14:15))];
% %
% % compute the total number of points per element
% npel = prod(lr1);
% %
% % compute number of active dimensions
% ndim = 2 + (lr1(3)>1);
% %
% % number of elements
% nel = str2double(header(17:26));
% %
% % number of elements in the file
% nelt = str2double(header(28:37));
% %
% % time
% time = str2double(header(39:58));
% %
% % istep
% istep = str2double(header(60:68));
% %
% % get file id  
% fid = str2double(header(70:75)); % TODO: multiple files not supported
% %
% % get tot number of files
% nf = str2double(header(77:82)); % TODO: multiple files not supported
% %
% % getfields [XUPTS]
% fields = strtrim(header(84:end));
% var=zeros(1,5);
% if sum(fields == 'X') > 0
%   var(1) = ndim;
% end
% if sum(fields == 'U') > 0
%   var(2) = ndim;
% end
% if sum(fields == 'P') > 0
%   var(3) = 1;
% end
% if sum(fields == 'T') > 0
%   var(4) = 1;
% end
% if sum(fields == 'S') > 0
%       var(5) = str2double(fields(end-1:end)); % TODO: scalars not implemented
% end
% nfields = sum(var);

