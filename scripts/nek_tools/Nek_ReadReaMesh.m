function [EL] = Nek_ReadReaMesh(fid,nelg,ndim,ifheat,npscal,nekver)

% INPUT: already open file handle for the rea file

disp('Reading mesh from rea file')

EL = [];

nsides=2^ndim;
Group = zeros(nelg);
for i=1:nelg
  tline = fgetl(fid);         % ELEMENT          8 [    1 ]    GROUP     0
  A =textscan(tline, '%s%d%s%d%s%s%d');   % Read element header
  GlElNo(i)=A{2};
  GroupNo(i)=A{7};

  tline = fgetl(fid);         % X1 X2 X3 X4
  xt = textscan(tline, '%f');   % Read X
  xtmp = xt{1};
  tline = fgetl(fid);         % Y1 Y2 Y3 Y4
  yt = textscan(tline, '%f');   % Read Y
  ytmp = yt{1};

  if (ndim==3)

    tline = fgetl(fid);         % Z1 Z2 Z3 Z4
    zt = textscan(tline, '%f');   % Read Z
    ztmp = zt{1};

    tline = fgetl(fid);         % X5 X6 X7 X8
    xt2 = textscan(tline, '%f');   % Read X
    tline = fgetl(fid);         % Y5 Y6 Y7 Y8
    yt2 = textscan(tline, '%f');   % Read Y
    tline = fgetl(fid);         % Z5 Z6 Z7 Z8
    zt2 = textscan(tline, '%f');   % Read Z
   
    xtmp = [xtmp; xt{1}];
    ytmp = [ytmp; yt{1}];
    ztmp = [ztmp; zt{1}];
    
    ZC(:,i) = ztmp;
  end
  XC(:,i) = xtmp;
  YC(:,i) = ytmp;

end

EL.Nelg = nelg;
EL.Ndim = ndim;
EL.GLobalNo = GlElNo;
EL.GroupNo  = GroupNo;
EL.XC = XC;
EL.YC = YC;
if (ndim==3)
  EL.ZC = ZC;
end  


%% 
disp('Reading Curved side data')
tline = fgetl(fid);  %***** CURVED SIDE DATA *****  
tline = fgetl(fid);  % 36 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE
cell =textscan(tline, '%f');
NCurve = cell{1};

for i=1:NCurve
  tline = fgetl(fid);             % Read curved data
  if (nelg<1000)
    iedge(i) = str2double(tline(1:3));
    ieg(i)   = str2double(tline(4:6));
    cell = textscan(tline(7:end),'%f %f %f %f %f %s');
  elseif (nelg<10^6)
    iedge(i) = str2double(tline(1:2));
    ieg(i)   = str2double(tline(3:8));
    cell = textscan(tline(9:end), '%f %f %f %f %f %s');
  else  
    iedge(i) = str2double(tline(1:2));
    ieg(i)   = str2double(tline(3:14));
    cell = textscan(tline(15:end),'%f %f %f %f %f %s');
  end

  % 5 parameters for describing the curve. Hard-coded
  for j=1:5
    cpar(j,1)=cell{j};
  end
  CurveParams(:,i) = cpar;
  CurveType(i) = cell{6};

end

EL.NCurve = NCurve;
EL.CurveEdeges = iedge;
EL.CurveIEG  = ieg;
EL.CurveParams = CurveParams;
EL.CurveType = CurveType;

%% Boundary Conditions
disp('Reading Boundary Conditions')
tline = fgetl(fid);    % ***** BOUNDARY CONDITIONS ***** 

nfldt = 1;
if (ifheat)
  nfldt = 2+npscal;
end
nbcs = nfldt;

for i=1:nbcs
  tline = fgetl(fid);    % ***** FLUID   BOUNDARY CONDITIONS *****
  cell = textscan(tline,'%s');
  tmp=max(strcmpi(cell{1},'NO'));
  ACTIVE=~tmp;

  if (ACTIVE)
    if nekver<=2.52
      nbcrea=3;
    elseif nekver>=2.55
      nbcrea=5;
    end
    rdfmt = [];
    for j=1:nbcrea
      rdfmt = [rdfmt '%f '];
    end

    ieg = [];

    for ieg= 1:nelg
      for iside=1:nsides
        tline = fgetl(fid);             % Read curved data
        if (nelg<1000)
          chtmp = tline(1);
          cbc = tline(2:4);
          I1 = str2double(tline(5:7));
          I2 = str2double(tline(8:10));
          cell=textscan(tline(11:end),rdfmt);
        elseif (nelg<10^5)
          chtmp = tline(1);
          cbc = tline(2:4);                        % Face BC
          I1 = str2double(tline(5:9));             % GLELno
          I2 = str2double(tline(10:10));           % side
          cell=textscan(tline(11:end),rdfmt);
        elseif (nelg<10^6)
          chtmp = tline(1);
          cbc = tline(2:4);
          I1 = str2double(tline(5:10));
%         I'm guessing its a bug that there is no I2.
%         That would mean the face number is missing.
%         So I'm adding it in
          I2 = str2double(tline(11:11));           % side
          cell=textscan(tline(12:end),rdfmt);
        else  
          chtmp = tline(1);
          cbc = tline(2:4);
          I1 = str2double(tline(5:12));
%         I'm guessing its a bug that there is no I2.
%         That would mean the face number is missing.
%         So I'm adding it in
          I2 = str2double(tline(13:13))           % side
          cell=textscan(tline(14:end),rdfmt);
        end   % nelg<

        CBC(I2,I1).BC = cbc;
        CBC(I2,I1).ConnectsTo = cell{1};
        CBC(I2,I1).OnFace = cell{2};
        CBC(I2,I1).PARAM3 = cell{3};
        CBC(I2,I1).PARAM4 = cell{4};
        CBC(I2,I1).PARAM5 = cell{5};

      end     % iside=1:nfaces
    end       % ieg=1:nelg    

  else      % Not active
    continue
  end       % if (ACTIVE)            

end         %i=1:nbcs 
if (nbcs==1)
  tline = fgetl(fid);    % ***** NO THERMAL BOUNDARY CONDITIONS *****
end  

EL.CBC = CBC;

%% Restart conditions
tline = fgetl(fid);  %  1 PRESOLVE/RESTART OPTIONS  *****
cell = textscan(tline, '%f');
Nrestart = cell{1};
for i=1:Nrestart
  tline = fgetl(fid);  % Restart file(s)
  cell=textscan(tline,'%s %s');
  rstFiles{i}=cell{1}{1};
  rstOptions{i}=cell{2}{1};
end

EL.rstFiles = rstFiles;
EL.rstOptions = rstOptions;




return

% Connectivity information in Nek
% <TYPE>  <Description>                      <Needed Parameters>        <Number of needed Parameters>
% E       internal (element connectivity)    adjacent element and face  2
% P       periodic                           periodic element and face  2 
% T       Dirichlet temperature/scalar       value                      1
% V       Dirichlet velocity                 u,v,w                      3
% O       outflow                            -                          0
% W       wall (no slip)                     -                          0                               
% F       flux                               flux                       1
% SYM     symmetry                           -                          0
% A       axisymmetric boundary              -                          0
% MS      moving boundary                    -                          0
% I       insulated (zero flux) for temperature                         0
% ON      Outflow, Normal (need surface to be normal to x,y,or z        0 


% Stored under Boundary Conditions in .rea file
%Examples:
%P      1       1      6.00000       3.00000       0.00000       0.00000       0.00000
%Face 1 of element 1 is periodic with element 6, face 3
%
%T       1       3      200.0         0.00000       0.00000       0.00000       0.00000
%Face 3 of element 1 has a constant value equal to 200.0
%
%t      1       3      0.00000       0.00000       0.00000       0.00000       0.00000
%Face 3 of element 1 will be assigned a value in USERBC 
