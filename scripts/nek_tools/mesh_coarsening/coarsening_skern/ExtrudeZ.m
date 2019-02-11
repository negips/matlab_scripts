function [XC,YC,ZC,GL,LayerGEl,EF] = ExtrudeZ(mesh2d,LayerGEl,glno,il,nz0,Lz,ifcl,cz_pl,zf1,zf2);

  il1   = il;  
  ind1  = mesh2d.layerindex{il1}; 
  LE1   = mesh2d.globalno(ind1); 

  nel_lay = length(ind1);

  XC  = []; YC  = []; ZC  = []; GL  = []; EF  = [];

  j=1;
% Go through all the elements in the 2D layer 
  while j<=nel_lay

    e=LE1(j);
    nz=nz0;
    dz=Lz/nz;

    etype = mesh2d.EType{e};
    elf1 = mesh2d.cbc(1,e).connectsto;    % Element no connecting face 1
    elf2 = mesh2d.cbc(2,e).connectsto;    % Element no connecting face 2
    elf3 = mesh2d.cbc(3,e).connectsto;    % Element no connecting face 3
    elf4 = mesh2d.cbc(4,e).connectsto;    % Element no connecting face 4
%   Element types of neighbors
    if elf1~=0
      f1t  = mesh2d.EType{elf1};
    else
      f1t  = 'B';                   % boundary
    end
    if elf2~=0
      f2t  = mesh2d.EType{elf2};
    else
      f2t  = 'B';
    end
    if elf3~=0
      f3t  = mesh2d.EType{elf3};
    else
      f3t  = 'B';
    end
    if elf4~=0
      f4t  = mesh2d.EType{elf4};
    else
      f4t  = 'B';
    end  

    [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_0(e,j,il,nz,dz,Lz,mesh2d,zf1,zf2);

    GL1  = GL1 + glno;
    
    XC = [XC XC1];
    YC = [YC YC1];
    ZC = [ZC ZC1];
    GL = [GL GL1];
    EF = [EF EF1];

    glno = max(GL);
    LayerGEl{e}=GL1;

    j=j+1;
      
  end


end   % function

 

