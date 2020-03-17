%  Script to read xml color table and create colormap

clear
clc

fname = 'blue-orange-div';
fname = 'gray-gold';
fname = 'mellow-rainbow';
fname = '5-step-melow-wave';
fname = '3-wave-yellow-grey-blue';
fname = '3w_bgYr';
fname = '4Wmed8';
fname = 'turqoise-olive';
fname = 'asym-orange-blue';
fname = 'green-brown-div';

xmlfile = [fname '.xml'];
xml = xml2struct(xmlfile);

np = length(xml.ColorMaps.ColorMap.Point);      % No of divisions
Pt = xml.ColorMaps.ColorMap.Point;

ctfile = ['00_' fname '.ct'];

fid = fopen(ctfile,'w');

node0 = '<?xml version="1.0"?>';
fprintf(fid,'%s\n',node0);

node1 = '<Object name="ColorTable">';
fprintf(fid,'%s\n',node1);

node2 = '<Field name="Version" type="string">2.12.0</Field>';
fprintf(fid,'\t%s\n',node2);

node3 = '<Object name="ColorControlPointList">';
fprintf(fid,'\t%s\n',node3);

% Scaling
sc = 255;

ll = 1:np;
if np>40
  ll=1:4:np;
elseif np>20
  ll=1:2:np;
end  

for i=ll
  nodei = '<Object name="ColorControlPoint">';
  fprintf(fid,'\t\t%s\n',nodei);

  r=round(str2num(Pt{i}.Attributes.r)*sc);
  g=round(str2num(Pt{i}.Attributes.g)*sc);
  b=round(str2num(Pt{i}.Attributes.b)*sc);
  fld1 = ['<Field name="colors" type="unsignedCharArray" length="4">',num2str(r,16),' ',num2str(g,16),' ', num2str(b,16),' ', num2str(1.0*sc),' </Field>'];
  fprintf(fid,'\t\t\t%s\n',fld1);

  fld2 = ['<Field name="position" type="float">',num2str((i-1)/(np-1)),'</Field>'];
  fprintf(fid,'\t\t\t%s\n',fld2);

  nodei = '</Object>';
  fprintf(fid,'\t\t%s\n',nodei);
end

node3='<Field name="category" type="string">UserDefined</Field>';
fprintf(fid,'\t\t%s\n',node3);

node2='</Object>';
fprintf(fid,'\t%s\n',node2);

node1='</Object>';
fprintf(fid,'%s\n',node1);

fclose(fid)

system(['mv ',ctfile, ' ~/.visit/' ])



