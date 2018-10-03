% Creating new structures/arrays so that future manipulation is easy

%clear
%clc
%close all
%
%load saab_wing2d.mat

% Create new Structures
% Simpler structure otherwise we end up with a very complicated code

for i=1:nlayers

  thislayer = LayersEl{i};
  thisfaceo = LayersFopO{i};
  thisfacev = LayersFopV{i};

  LX = [];
  LY = [];
  LE = [];
  LBC= [];

  nel = length(thislayer);
  for j=1:nel

%   First layer first element is complicated.
    e=thislayer(j);
    fo=thisfaceo(j);
    fv=thisfacev(j);

    [xcs ycs bcs] = GetFirstEl(rea,e,fo,fv);
    LX = [LX xcs];
    LY = [LY ycs];
    LE = [LE e];
    LBC= [LBC bcs];

  end

  LayerX{i}=LX;
  LayerY{i}=LY;
  LayerE{i}=LE;
  LayerBC{i}=LBC;
 
end


function [xcs ycs bcs] = GetFirstEl(rea,e,fo,fv)

    dtol=1.0e-12;     % Distance tolerance
      
    fop = fo;  
    bc3 = rea.mesh.cbc(fop,e).bc;  

    iop = [fop fop+1];
    if (iop(2)>4)
      iop(2)=1;
    end

    fo=fo+2;
    if fo>4
      fo=fo-4;
    end
    bc1 = rea.mesh.cbc(fo,e).bc;  

    io = [fo fo+1];
    if (io(2)>4)
      io(2)=1;
    end  

    fov=fv;  
    bc2 = rea.mesh.cbc(fov,e).bc;  

    fv=fv+2;
    if fv>4
      fv=fv-4;
    end  
    bc4 = rea.mesh.cbc(fv,e).bc;  

    iv = [fv fv+1];
    if (iv(2)>4)
      iv(2)=1;
    end  
   
    xo = rea.mesh.xc(io,e);
    yo = rea.mesh.yc(io,e);

    xv = rea.mesh.xc(iv,e);
    yv = rea.mesh.yc(iv,e);

    d1 = sqrt((xo-xv).^2 + (yo-yv).^2);
    d2 = sqrt((xo-flipud(xv)).^2 + (yo-flipud(yv)).^2);

    if abs(d1(1))< dtol || abs(d2(1))< dtol
      x1 = xo(1);
      y1 = yo(1);
      x2 = xo(2);
      y2 = yo(2);

      if abs(d1(1))< dtol
        x4 = xv(2);
        y4 = yv(2);
      else  
        x4 = xv(1);
        y4 = yv(1);
      end

    else  
      x1 = xo(2);
      y1 = yo(2);
      x2 = xo(1);
      y2 = yo(1);

      if abs(d1(2))< dtol
        x4 = xv(1);
        y4 = yv(1);
      else  
        x4 = xv(2);
        y4 = yv(2);
      end


    end  

    xop = rea.mesh.xc(iop,e);
    yop = rea.mesh.yc(iop,e);

    d3 = sqrt((xop-x4).^2 + (yop-y4).^2);
    if abs(d3(1)) < dtol
      x3 = xop(2);
      y3 = yop(2);
    else  
      x3 = xop(1);
      y3 = yop(1);
    end

    xcs = [x1; x2; x3; x4];
    ycs = [y1; y2; y3; y4];
    bcs = [bc1; bc2; bc3; bc4]; 


end


