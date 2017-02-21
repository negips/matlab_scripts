%% read all experimental data files from folder 


[status,result] = system(['ls ' folder '*']);

inds1 = strfind(result,fol);
inds2 = strfind(result,'.h5');

nfiles = length(inds2);
U0=zeros(nfiles,1);
alpha=zeros(nfiles,1)-99;
defl=zeros(nfiles,1)-99;
ifturb=zeros(nfiles,1);

disp(['N files: ', num2str(nfiles)])

for i=1:nfiles

    ind1 = inds1(i)+lfol;
    ind2 = inds2(i)+2;
    fname=result(ind1:ind2);
    filenames{i}=fname;

    inds4=strfind(fname,'alphasweep');
    inds3=strfind(fname,'_');

    if ~isempty(inds4)
      ind1=2;
      ind2=inds3(1)-1;      
      U0(i) = str2double(fname(ind1:ind2));
      ind1=inds3(1)+2;
      ind2=inds3(2)-1;
      defl(i) = str2double(fname(ind1:ind2));
      if isnan(defl(i)) && i>1 
        defl(i)=defl(i-1);
      end  
      continue
    end
   
    if (length(inds3)>2)
      % Freestream
      ind1=2;
      ind2=inds3(1)-1;
      U0(i) = str2double(fname(ind1:ind2));
      % apha
      ind1=inds3(1)+2;
      ind2=inds3(2)-1;
      alpha(i) = str2double(fname(ind1:ind2));
      % flap deflection      
      ind1=inds3(2)+2;
      ind2=inds3(3)-1;
      defl(i) = str2double(fname(ind1:ind2));
    else
      % Freestream
      ind1=2;
      ind2=inds3(1)-1;
      U0(i) = str2double(fname(ind1:ind2));
      ind1=inds3(1)+1;
      ind2=inds3(2)-1;
      % flap deflection
      defl(i) = str2double(fname(ind1:ind2));
    end

    inds4=strfind(fname,'turb');
    if ~isempty(inds4)
      ifturb(i)=1;
    end  
    
end

