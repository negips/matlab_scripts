function [XC,YC,ZC,GL] = QuadExpansion(mesh2d,il,nz0,Lz);

  il1   = il;  
  ind1  = mesh2d.layerindex{il1}; 
  LE1   = mesh2d.globalno(ind1); 

  nel_lay = length(ind1);

  XC  = []; YC  = []; ZC  = []; GL  = [];

  j=1;
  while j<=nel_lay

    XC1 =[]; YC1 =[]; ZC1 =[]; GL1 =[];
    XC2 =[]; YC2 =[]; ZC2 =[]; GL2 =[];
    XC3 =[]; YC3 =[]; ZC3 =[]; GL3 =[];
    XC4 =[]; YC4 =[]; ZC4 =[]; GL4 =[];

    e=LE1(j);
    nz=nz0;

    etype = mesh2d.EType{e};
    elf1 = mesh2d.cbc(1,e).connectsto;    % Element no connecting face 1
    elf3 = mesh2d.cbc(3,e).connectsto;    % Element no connecting face 3
    elf4 = mesh2d.cbc(4,e).connectsto;    % Element no connecting face 4
%   Element types of neighbors
    if elf1~=0
      f1t  = mesh2d.EType{elf1};
    else
      f1t  = 'B';                   % boundary
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

      [XC1,XC2,YC1,YC2,ZC1,ZC2,GL1,GL2]=BuildLayer_S(XC1,XC2,YC1,YC2,ZC1,ZC2,GL1,GL2,e,j,il,nz,Lz,mesh2d);

    elseif strcmpi(etype,'e3') && strcmpi(f3t,'e4')
%     done in first if condition
    elseif strcmpi(etype,'e1') && strcmpi(f1t,'e4')
%     done in first if condition

    elseif strcmpi(etype,'e4') && strcmpi(f4t,'e4')       % Left elongated element 

      [XC1,YC1,ZC1,GL1,j]=BuildLayer_EL(XC1,YC1,ZC1,GL1,e,j,il,nz,Lz,mesh2d);

    elseif strcmpi(etype,'e4') && strcmpi(f4t,'e1')       % Right elongated element 

      [XC1,YC1,ZC1,GL1,j]=BuildLayer_ER(XC1,YC1,ZC1,GL1,e,j,il,nz,Lz,mesh2d);

    else  
       disp(['Some more elements need attention ', etype,' ', f1t,' ', f3t,' ', f4t])
       figure(4);
       plot(mesh2d.xc(:,e),mesh2d.yc(:,e)); hold on
    end

    j=j+1;

    XC = [XC XC1 XC2 XC3 XC4];
    YC = [YC YC1 YC2 YC3 YC4];
    ZC = [ZC ZC1 ZC2 ZC3 ZC4];
    GL = [GL GL1 GL2 GL3 GL4];

  end             % j<=nel_layer
      
end         % function
%---------------------------------------------------------------------- 

