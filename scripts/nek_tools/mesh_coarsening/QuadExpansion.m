function [XC,YC,ZC,GL,LayerGEl,EF] = QuadExpansion(mesh2d,LayerGEl,glno,il,nz0,Lz,ifcl,cz_pl,zf1,zf2);

  il1   = il;  
  ind1  = mesh2d.layerindex{il1}; 
  LE1   = mesh2d.globalno(ind1); 

  nel_lay = length(ind1);

  XC  = []; YC  = []; ZC  = []; GL  = []; EF  = [];

  XC1 = []; YC1 = []; ZC1 = []; GL1 = []; EF1 = [];

  j=1;
% Go through all the elements in the 2D layer 
  while j<=nel_lay

    e=LE1(j);
    nz=nz0;

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

    if strcmpi(etype,'s') 
      [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_S(e,j,il,nz,Lz,mesh2d,zf1,zf2);

%    elseif (strcmpi(etype,'e1') || strcmpi(etype,'e3'))
%      [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_E13(e,j,il,nz,Lz,mesh2d,zf1,zf2);

    elseif strcmpi(etype,'e4') && strcmpi(f3t,'e4')       % Left elongated element 

      [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_EL(e,j,il,nz,Lz,mesh2d,zf1,zf2);

    elseif strcmpi(etype,'e4') && strcmpi(f1t,'e4')       % Right elongated element 

      [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_ER(e,j,il,nz,Lz,mesh2d,zf1,zf2);

    elseif strcmpi(etype,'e4') && strcmpi(f3t,'B')        % Elongated Boundary element 

      [XC1,YC1,ZC1,GL1,EF1]=BuildLayer_EL(e,j,il,nz,Lz,mesh2d,zf1,zf2);

    else  
       disp(['Some more elements need attention ', etype,' ', f1t,' ', f2t, ' ', f3t,' ', f4t])
       figure(4);
       xt=[mesh2d.xc(:,e); mesh2d.xc(1,e)]; 
       yt=[mesh2d.yc(:,e); mesh2d.yc(1,e)]; 

       plot(xt,yt); hold on
    end

    j=j+1;

    GL1 = GL1+glno;

    XC = [XC XC1];
    YC = [YC YC1];
    ZC = [ZC ZC1];
    GL = [GL GL1];
    EF = [EF EF1];

    glno=max(GL);

    LayerGEl{e}=GL1;


  end             % j<=nel_layer
      
end         % function
%---------------------------------------------------------------------- 

