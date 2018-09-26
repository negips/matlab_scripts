function [REA] = Nek_ReadRea(casename)

% Function to read rea files (Mostly for the mesh right now)
% Maybe I should directly write a reader for re2

casename = 'lu';
reafile = [casename '.rea'];

fid=fopen(reafile,'r');

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
  tline = fgetl(fid); 
  cell  = sscanf(tline(2:end), '%s');
  Logical{i,1} = cell(1);
  Logical{i,2} = cell(2:end);

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
  EL = re2_mesh_read(casename,Nelgv,Ndim,IFHEAT,NPSCAL);
else
  EL = Nek_ReadReaMesh(fid,Nelgv,Ndim,IFHEAT,NPSCAL,NekVer);
end  

EL.xfac  = Xfac;
EL.yfac  = Yfac;
EL.xzero = Xzero;
EL.yzero = Yzero;
EL.ifre2 = IFRE2;
EL.ifgtp = IFGTP;

REA.casename = casename;
REA.nekver   = NekVer;
REA.nparams  = nparams;
REA.nlogical = nLogic;
REA.npscal   = NPSCAL;
REA.ifflow   = IFFLOW;
REA.ifheat   = IFHEAT;
REA.param    = PARAM;
REA.logical  = Logical;
REA.mesh     = EL;

fclose(fid);

return
