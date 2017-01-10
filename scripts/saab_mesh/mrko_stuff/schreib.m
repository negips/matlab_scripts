function status = schreib(data,nel,header,tag,cnttail,dim,endian,fname,fields)

%%% This function writes binary data in the nek5000 new file format

% Open file
[fid,message] = fopen(fname,'w+',['ieee-' endian]);

if fid == -1
    disp(message) 
    return 
end

% New header
header(20:39) = '                    ';
snel = num2str(nel);
header(23:22+length(snel)) = snel;
header(34:33+length(snel)) = snel;

fwrite(fid,header,'*char');

fwrite(fid,tag,'*float64');

% Element numbers
lglel = 1:nel;
lglel = lglel'; 
fwrite(fid,lglel,'*int32');


% How many fields?
if strcmp(fields,'xup')
    if dim == 2
        nfields = 6;
    elseif dim == 3
        nfields = 8;
    end
elseif strcmp(fields,'up')
    if dim == 2
        nfields = 3;
    elseif dim == 3
        nfields = 4;
    end
end


% Write data
for iel=1:nel
    for ifld = 1:dim
        fwrite(fid,data(iel,:,ifld),'*float64'); %1:N^dim
    end
end
for iel=1:nel
    for ifld = dim+1:2*dim
        fwrite(fid,data(iel,:,ifld),'*float64');
    end
end
for iel=1:nel
    fwrite(fid,data(iel,:,nfields-1),'*float64');
end
for iel=1:nel
    fwrite(fid,data(iel,:,nfields),'*float64');
end
% % Not understood: More data to be written in the fld file.
% % Pad with zeros.
% cnttail = nel*cnttail;
% quatsch = zeros(1,cnttail);
% fwrite(fid,quatsch,'*float32');

status = fclose(fid);

return
