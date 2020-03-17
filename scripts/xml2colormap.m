%  Script to read xml color table and create colormap

clear
clc

i=1;
fname{i} = 'blue-orange-div';
i=i+1;
fname{i} = 'gray-gold';
i=i+1;
fname{i} = 'mellow-rainbow';
i=i+1;
fname{i} = '5-step-melow-wave';
i=i+1;
fname{i} = '3-wave-yellow-grey-blue';
i=i+1;
fname{i} = '3w_bgYr';
i=i+1;
fname{i} = '4Wmed8';
i=i+1;
fname{i} = 'turqoise-olive';
i=i+1;
fname{i} = 'asym-orange-blue';
i=i+1;
fname{i} = 'green-brown-div';


nfiles=length(fname);

for i=1:nfiles
  file=fname{i};
  xmlfile = [file '.xml'];
  xml = xml2struct(xmlfile);
  
  np = length(xml.ColorMaps.ColorMap.Point);      % No of divisions
  Pt = xml.ColorMaps.ColorMap.Point;
  
  matfile = ['00_' file '.mat'];
  
  % Scaling
  sc = 1.0;
  ll = 1:np;
  
  colormap=[];
  for i=ll
    r=str2num(Pt{i}.Attributes.r)*sc;
    g=str2num(Pt{i}.Attributes.g)*sc;
    b=str2num(Pt{i}.Attributes.b)*sc;
  
    colormap = [colormap; r g b];
  end
  save(matfile, 'colormap')

end  


