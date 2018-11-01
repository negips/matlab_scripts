function [REA] = Nek_ReadRea(casename)

% Function to read rea files (Mostly for the mesh right now)
% Maybe I should directly write a reader for re2

reafile = [casename '.rea'];

fid=fopen(reafile,'r');

if fid == -1
  disp(message) 
  return 
end

disp(['Reading ' reafile])

tline = fgetl(fid);                                              % ***** Parameters ****
tline = fgetl(fid); cell=textscan(tline, '%f');                  % Nekton Version
NekVer = cell{1};
tline = fgetl(fid); cell=textscan(tline, '%f');                  % 2D/3D Simulation
Ndim = cell{1};
tline = fgetl(fid); cell=textscan(tline, '%f');                  % No of parameters
nparams = cell{1};

% Should probably create a cell array for all parameter descriptions
PARAM = zeros(200,1);
for i=1:nparams
  tline = fgetl(fid); cell=textscan(tline, '%f');                % Param value
  PARAM(i) = cell{1};
end
NPSCAL = round(PARAM(23));                                       % No. of Passive scalars

% npsLines of passive scalar data
tline = fgetl(fid); cell=textscan(tline, '%f');                  % Passive scalar data
npsLines = cell{1};
if (NPSCAL>0 && npsLines>0)
  % Skipping this for now...  
else
  % Skip these lines
  for i=1:npsLines
    tline = fgetl(fid); cell=textscan(tline, '%f');              % Passive scalar lines
    % Not saving into anything right now
  end
end   

tline = fgetl(fid); cell=textscan(tline, '%f');                  % Logical Switches
nLogic = cell{1};
for i=1:nLogic
  tline  = fgetl(fid); 
  cell   = sscanf(tline(2:end), '%s');
  ifnav  = strfind(lower(cell),'ifnav');
  ifadvc = strfind(lower(cell),'ifadvc');
  iftmsh = strfind(lower(cell),'iftmsh');

  if ~isempty(ifnav) && ~isempty(ifadvc)
%   Navier Stokes and passive scalar advection switches
    indmin = min([ifnav ifadvc]);
    Logical{i,1} = cell(1:indmin-1);
    Logical{i,2} = 'IFNAV && IFADVC';
  elseif ~isempty(ifnav) && isempty(ifadvc) 
    indmin = ifnav;
    Logical{i,1} = cell(1:indmin-1);
    Logical{i,2} = 'IFNAV';
  elseif isempty(ifnav) && ~isempty(ifadvc) 
    indmin = ifadvc;
    Logical{i,1} = cell(1:indmin-1);
    Logical{i,2} = 'IFADVC';
  elseif ~isempty(iftmsh) 
    indmin = iftmsh;
    Logical{i,1} = cell(1:indmin-1);
    Logical{i,2} = 'IFTMSH';
  else
    Logical{i,1} = cell(1);
    Logical{i,2} = cell(2:end);
  end  

% Checking for IFFLOW and IFHEAT logicals.
% IfHEAT is required for reading the mesh
  tmp1 = strfind(Logical{i,2},'ifflow');
  tmp2 = strfind(Logical{i,2},'IFFLOW');
  if ~isempty(tmp1) || ~isempty(tmp2)
    if (strcmpi(Logical{i,1},'T'))
      IFFLOW = true;
    else  
      IFFLOW = false;
    end
  end  

% IFHEAT  
  tmp1 = strfind(Logical{i,2},'ifheat');
  tmp2 = strfind(Logical{i,2},'IFHEAT');
  if ~isempty(tmp1) || ~isempty(tmp2)
    if (strcmpi(Logical{i,1},'T'))
      IFHEAT = true;
    else  
      IFHEAT = false;
    end
  end  

end % i=1:nLogic

tline = fgetl(fid); cell=textscan(tline, '%f');                  % Mesh parameter
% These are basicly for pretex (I think)
Xfac=cell{1}(1); 
Yfac=cell{1}(2); 
Xzero=cell{1}(3); 
Yzero=cell{1}(4);

tline = fgetl(fid);                                              % **Mesh Data**
tline = fgetl(fid); cell=textscan(tline, '%f');                  % Mesh parameter
Nelgs = cell{1}(1);
meshndim  = cell{1}(2);
Nelgv = cell{1}(3);

IFRE2 = false;                                                       % If .re2 file
if (Nelgs<0) 
  IFRE2=true;
end

IFGTP =false;                                                       % if Global tensor product
if (meshndim<0)
  IFGTP = true;
end  

if (IFRE2)
%  Not implemented      
%  EL = re2_mesh_read(casename,Nelgv,Ndim,IFHEAT,NPSCAL);
else
  EL = Nek_ReadReaMesh(fid,Nelgv,Ndim,IFHEAT,NPSCAL,NekVer);
end  

EL.xfac  = Xfac;
EL.yfac  = Yfac;
EL.xzero = Xzero;
EL.yzero = Yzero;
EL.ifre2 = IFRE2;
EL.ifgtp = IFGTP;

%% Restart conditions
tline = fgetl(fid);  %  1 PRESOLVE/RESTART OPTIONS  *****
cell = textscan(tline, '%f');
Nrestart = cell{1};
rstFiles=[];
rstOptions=[];
for i=1:Nrestart
  tline = fgetl(fid);  % Restart file(s)
  cell=textscan(tline,'%s %s');
  rstFiles{i}=cell{1}{1};
  rstOptions{i}=cell{2}{1};
end

tline = fgetl(fid); cell=textscan(tline, '%f');                % Initial conditions
nic=cell{1};
initialconditions=[];
for i=1:nic
  initialconditions{i} = fgetl(fid);
end  

tline = fgetl(fid);                                            % ****Drive Force Data****
tline = fgetl(fid); cell=textscan(tline, '%f');                % n lines
ndriveforce = cell{1};
driveforce=[];
for i=1:ndriveforce
  driveforce{i} = fgetl(fid);
end

% This may not be entirely correct. Needs to be checked
tline = fgetl(fid);                                            % ****Variable Property Data****
tline = fgetl(fid); cell=textscan(tline, '%f');                % n lines
nvplines = cell{1};

tline = fgetl(fid); cell=textscan(tline, '%f');                % n Packets of data
npackets  = cell{1};
datapacket=[];
for i=1:npackets
  datapacket{i}=  fgetl(fid);                                  % Not sure what this is
end
%

tline = fgetl(fid);                                            % ****History and Integral Data****
tline = fgetl(fid); cell=textscan(tline, '%f');                % n points
nhist = cell{1};
history=[];
for i=1:nhist
  history{i}=  fgetl(fid);
end

tline = fgetl(fid);                                            % ****Output field specification****
tline = fgetl(fid); cell=textscan(tline, '%f');                % n specifications
noutspec = cell{1};
outputspec=[];
for i=1:noutspec
  tline=  fgetl(fid);
  kps=strfind(lower(tline), 'passive scalars');
  ktg=strfind(lower(tline), 'temperature gradient');
  if ~isempty(kps)
    cell = textscan(tline,'%f');
    outputspec{i,1}=num2str(cell{1});
    index=strfind(tline,num2str(cell{1}));
    str  =strtrim(tline(index+1:end));
    outputspec{i,2}=str;
  elseif ~isempty(ktg)
    cell = textscan(tline,'%s%s%s');
    outputspec{i,1}=cell{1}{1};
    outputspec{i,2}=[cell{2}{1} ' ' cell{3}{1}];
  else
    cell   = sscanf(tline, '%s');
    outputspec{i,1}=cell(1);
    outputspec{i,2}=cell(2:end);
  end
end

tline = fgetl(fid);                                            % ****Object specification****
objects=[];
i=0;
while ~feof(fid)
  tline=  fgetl(fid);
  i=i+1;
  cell  = sscanf(tline, '%f');
  objects{i,1}=cell;
  index=strfind(tline,num2str(cell));
  str=strtrim(tline(index+1:end));   % Remove leading and trailing white spaces
  objects{i,2}=str;

end
nobjects=i;

REA.casename = casename;
REA.nekver   = NekVer;
REA.nparams  = nparams;
REA.Nlogical = nLogic;
REA.npscal   = NPSCAL;
REA.ifflow   = IFFLOW;
REA.ifheat   = IFHEAT;
REA.param    = PARAM;
REA.logical  = Logical;
REA.mesh     = EL;
REA.Nrestart = Nrestart;
REA.rstFiles = rstFiles;
REA.rstOptions = rstOptions;
REA.Nic=nic;
REA.initialconditions  = initialconditions;
REA.Ndriveforce=ndriveforce;
REA.driveforce=driveforce;
REA.Nvplines=nvplines;
%REA.variableproperty=variableproperty;
REA.Npackets=npackets;
REA.datapacket=datapacket;
REA.Nhist=nhist;
REA.history=history;
REA.Noutspec=noutspec;
REA.outputspec=outputspec;
%REA.Npsdata=npsdata;
%REA.passivescalar=passivescalar;
REA.Nobjects=nobjects;
REA.objects=objects;


fclose(fid);

return
