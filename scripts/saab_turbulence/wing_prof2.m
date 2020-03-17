close all; 
clear all; 
clc
double precision;
format long;


for ssf=1:79
  tic
      
  streamwise = 1;
  svfname = sprintf('%s%4.4d%s', 'saab750k_pitch',ssf,'.mat');
  display(['Streamwise direction: ', num2str(streamwise)])
  display(['Save file name ' svfname])
  display('Calculating fields')
  
  % Fields in the binary record from stat files
  
  %Statistic fields
  % u,v,w,p are instantaneous quantities 
  % averaged in time and homogeneous direction z 
  
  % 1.  <u>           % F1
  % 2.  <v>           % F2
  % 3.  <w>           % F3
  % 4.  <p>           % F4
  
  % 5.  <uu>          % F5
  % 6.  <vv>          % F6
  % 7.  <ww>          % F7
  % 8.  <pp>          % F8
  
  % 9.  <uv>          % F9
  % 10. <vw>          % F10
  % 11. <uw>          % F11
  
  % 12. <pu>          % F12
  % 13. <pv>          % F13
  % 14. <pw>          % F14
  
  % 15. <pdudx>       % F15
  % 16. <pdudy>       % F16
  % 17. <pdudz>       % F17
  
  % 18. <pdvdx>       % F18
  % 19. <pdvdy>       % F19
  % 20. <pdvdz>       % F20
  
  % 21. <pdwdx>       % F21
  % 22. <pdwdy>       % F22
  % 23. <pdwdz>       % F23
  
  % 24. <uuu>         % F24
  % 25. <vvv>         % F25
  % 26. <www>         % F26
  % 27. <ppp>         % F27
  
  % 28. <uuv>         % F28
  % 29. <uuw>         % F29
  % 30. <vvu>         % F30
  % 31. <vvw>  	  % F31   
  % 32. <wwu>         % F32
  % 33. <wwv>         % F33
  % 34. <uvw>         % F34
  
  % 35. <uuuu>        % F35
  % 36. <vvvv>        % F36
  % 37. <wwww>        % F37
  % 38. <pppp>        % F38
  
  % 39. <uuuv>        % F39
  % 40. <uuvv>        % F40
  % 41. <uvvv> 	    % F41 
     
  % 42. e11: <((du/dx)^2+(du/dy)^2+(du/dz)^2)> % F42 
  % 43. e22: <((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)> % F43 
  % 44. e33: <((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)> % F44 
  % 45. e12: <(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F45  
  % 46. e13: <(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F46 
  % 47. e23: <(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F47 
  
  % 48. <dw/dx*dw/dx> % F48 
  % 49. <dw/dy*dw/dy> % F49 
  % 50. <dw/dx*dw/dy> % F50 
  
  % 51. <du/dx*du/dx> % F51
  % 52. <du/dy*du/dy> % F52 
  % 53. <du/dx*du/dy> % F53 
  
  % 54. <dv/dx*dv/dx> % F54
  % 55. <dv/dy*dv/dy> % F55 
  % 56. <dv/dx*dv/dy> % F56 
  
  % 57. <du/dx*dv/dx> % F57
  % 68. <du/dy*dv/dy> % F58 
  % 59. <du/dx*dv/dy> % F59 
  % 60. <du/dy*dv/dx> % F60
  
  %Derivative fields
  % 1. dU/dx           % D1
  % 2. dU/dy           % D2
  % 3. dV/dx           % D3
  % 4. dV/dy           % D4
  
  % 5. dW/dx           % D5
  % 6. dW/dy           % D6
  % 7. dP/dx           % D7
  % 8. dP/dy           % D8
  
  % 9.  d<uu>/dx       % D9
  % 10. d<uu>/dy       % D10
  % 11. d<vv>/dx       % D11
  % 12. d<vv>/dy       % D12
  
  % 13. d<ww>/dx       % D13
  % 14. d<ww>/dy       % D14
  % 15. d<pp>/dx       % D15
  % 16. d<pp>/dy       % D16
  
  % 17. d<uv>/dx       % D17
  % 18. d<uv>/dy       % D18
  % 19. d<vw>/dx       % D19
  % 20. d<vw>/dy       % D20
  
  % 21. d<uw>/dx       % D21
  % 22. d<uw>/dy       % D22
  % 23. d<uuu>/dx      % D23
  % 24. d<uuu>/dy      % D24
  
  % 25. d<vvv>/dx      % D25
  % 26. d<vvv>/dy      % D26
  % 27. d<www>/dx      % D27
  % 28. d<www>/dy      % D28
  
  % 29. d<ppp>/dx      % D29
  % 30. d<ppp>/dy      % D30
  % 31. d<uuv>/dx      % D31
  % 32. d<uuv>/dy      % D32
  
  % 33. d<uuw>/dx      % D33
  % 34. d<uuw>/dy      % D34
  % 35. d<vvu>/dx      % D35
  % 36. d<vvu>/dy      % D36
  
  % 37. d<vvw>/dx      % D37
  % 38. d<vvw>/dy      % D38
  % 39. d<wwu>/dx      % D39
  % 40. d<wwu>/dy      % D40
  
  % 41. d<wwv>/dx      % D41
  % 42. d<wwv>/dy      % D42
  % 43. d<uvw>/dx      % D43
  % 44. d<uvw>/dy      % D44
  
  % 45. d2U/dx2        % D45
  % 46. d2U/dy2        % D46
  % 47. d2V/dx2        % D47
  % 48. d2V/dy2        % D48
  
  % 49. d2W/dx2        % D49
  % 50. d2W/dy2        % D50
  % 51. d2<uu>/dx2     % D51
  % 52. d2<uu>/dy2     % D52
  
  % 53. d2<vv>/dx2     % D53
  % 54. d2<vv>/dy2     % D54
  % 55. d2<ww>/dx2     % D55
  % 56. d2<ww>/dy2     % D56
  
  % 57. d2<uv>/dx2     % D57
  % 58. d2<uv>/dy2     % D58
  % 59. d2<uw>/dx2     % D59
  % 60. d2<uw>/dy2     % D60
  
  % 61. d2<vw>/dx2     % D61
  % 62. d2<vw>/dy2     % D62
  
  % 63. d<pu>/dx       % D63
  % 64. d<pu>/dy       % D64
  % 65. d<pv>/dx       % D65
  % 66. d<pv>/dy       % D66
  
  % 67. d<pw>/dx       % D67
  % 68. d<pw>/dy       % D68
  
  
  %Generate geometry of NACA profile
  % [xu,yu,alphau,xl,yl,alphal]=naca_prof(c,x,t,m,p);
  
  % Read interpolated field: int_fld file
  %fname         = 'ZSTAT/int_fld004';
  fname = sprintf('%s%4.4d','ZSTAT/int_fld',ssf);
  [fid,message] = fopen(fname,'r','ieee-le');
  dum1          = fread(fid,1,'int32');
  F             = fread(fid,dum1,'*char')'  ;
  dum2          = fread(fid,2,'int32') ;
  
  Rer           = fread(fid,1,'*float64')   ;
  Domain        = fread(fid,3,'*float64')   ;
  mnelxyz       = fread(fid,3,'int32')      ;
  nel           = mnelxyz(1)                ;
  Poly          = fread(fid,3,'int32')      ;
  nstat         = fread(fid,1,'int32')      ;    
  nderiv        = fread(fid,1,'int32')      ;    
  times         = fread(fid,1,'*float64')   ;
  timee         = fread(fid,1,'*float64')   ;
  atime         = fread(fid,1,'*float64')   ;
  DT            = fread(fid,1,'*float64')   ;
  nrec          = fread(fid,1,'int32')      ;
  Tint          = fread(fid,1,'*float64')   ;
  npoints       = fread(fid,1,'int32')      ;
  wallpts       = fread(fid,1,'int32')      ;
  dum2          = fread(fid,2,'int32')   ;
  %dum6          = fread(fid,1,'*float64')   ;
  
  % Read interpolating mesh: x.fort and y.fort
  %fname_x           = 'ZSTAT/x.fort';
  %[fid_x,message_x] = fopen(fname_x,'r','ieee-le');
  %hdr_x             = fread(fid_x,1,'int32');
  %dum1              = fread(fid_x,3,'int32');
  %x_pts             = fread(fid_x,npoints,'*float64');
  %fclose(fid_x);
  %
  %fname_y          = 'ZSTAT/y.fort';
  %[fid_y,message_y] = fopen(fname_y,'r','ieee-le');
  %hdr_y             = fread(fid_y,1,'int32');
  %dum2              = fread(fid_y,3,'int32');
  %y_pts             = fread(fid_y,npoints,'*float64');
  %fclose(fid_y);
  
  %Find wall-normal vector and its length
  %in=find(x_pts==-1);
  %yn=y_pts(in);
  
  %Number of profiles on top and bottom
  np=wallpts;
  nypts=npoints/wallpts;  % wall normal points at each location
  ln=nypts;
  
  %
  xa              = fread(fid,wallpts,'*float64');
  dum1            = fread(fid,2,'int32');
  ya              = fread(fid,wallpts,'*float64');
  dum1            = fread(fid,2,'int32');
  
  snx             = fread(fid,wallpts,'*float64');
  dum1            = fread(fid,2,'int32');
  sny             = fread(fid,wallpts,'*float64');
  dum1            = fread(fid,2,'int32');
  
  %Define structures and coordinates
  top=struct;
  bottom=struct;
  
  tc=0;
  bc=0;
  side=-1*ones(wallpts,1);            % 1=> top. 0=> bottom
  for i=1:wallpts
    yslope=sny(i);
    if (yslope>=0)      % top
      side(i)=1;
      tc=tc+1;
      top(tc).xa=xa(i);
      top(tc).ya=ya(i);
      top(tc).snx=snx(i);
      top(tc).sny=sny(i);
      top(tc).alpha=atan2(sny(i),snx(i))-pi/2;          % y-axis is 0 degrees rotation    
    else
      side(i)=0;
      bc=bc+1;
      bottom(bc).xa=xa(i);
      bottom(bc).ya=ya(i);
      bottom(bc).snx=snx(i);
      bottom(bc).sny=sny(i);
      bottom(bc).alpha=atan2(sny(i),snx(i))-pi/2;
    end
  end
  topcount=tc;
  botcount=bc;
  tindicies=find(side);
  bindicies=find(~side); 
  
  % Read x values
  for i=1:wallpts
    xprof=fread(fid,nypts,'*float64');
    if side(i)==1
      pos=sum(side(1:i));
      top(pos).x=xprof;
    else
      pos=sum(~side(1:i));
      bottom(pos).x=xprof;
    end
  end          
  dum2=fread(fid,2,'int32');
  
  % Read y values
  for i=1:wallpts
    yprof=fread(fid,nypts,'*float64');
    if side(i)==1
      pos=sum(side(1:i));
      top(pos).y=yprof;
    else
      pos=sum(~side(1:i));
      bottom(pos).y=yprof;
    end
  end          
  dum2=fread(fid,2,'int32');
  
  
  % prabal
  % Assuming yn is the normal distance vector.
  for tc=1:topcount
    xv=top(tc).x;
    yv=top(tc).y;
    x0=xv(1);
    y0=yv(1);
    yn=sqrt((xv-x0).^2 + (yv-y0).^2);
    top(tc).yn=yn;
  end
  
  for bc=1:botcount
    xv=bottom(bc).x;
    yv=bottom(bc).y;
    x0=xv(1);
    y0=yv(1);
    yn=sqrt((xv-x0).^2 + (yv-y0).^2);
    bottom(bc).yn=yn;
  end
  
  %for i=1:np
  %    top(i).alpha=alphau(i);
  %    bottom(i).alpha=alphal(i);
  %    top(i).yn=yn;
  %    bottom(i).yn=yn;
  %end
  
  %Define data arrays
  Mat=zeros(ln,np);
  D=F;
  
  % Arrange stat fields in array form
  for ii=1:nstat+nderiv
      
      if ii == 1
          fseek(fid,0,'cof');
      else
          fseek(fid,8,'cof');
      end
      
      % Store current field as an array
      MatR = fread(fid,npoints,'*float64'); 
  
      % Arrange current field in array form
      for i=1:np
          Mat(:,i)=MatR(ln*(i-1)+1:ln*(i-1)+ln);
      end
          
      if ii<=nstat
          % Generata variable name for field X as FX
          v=genvarname('F', who);
      else
          % Generata variable name for field X as DX
          v=genvarname('D', who);
      end
      
      % Store matrix in variable FX or Dx
      evalc([v '=Mat']);   
      
  end
  
  fclose(fid);
  
  % Simulation parameters
  % Bulk Reynolds number Reb=Ub*c/nu
  Reb=Rer; 
  
  % Domain dimensions
  Lx=Domain(1);
  Ly=Domain(2);
  
  % Fluid density
  rho=1; 
  
  % Kinematic viscosity. Both Ub and c are unit normalizing parameters
  nu=1/Reb;
  
  % Molecular visosity
  mu=rho*nu;
  
  %Mean velocities. Tensors of Rank 1.
  tc=0;
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[F1(i,j); F2(i,j); F3(i,j)];
        top(tc).U(i,1)=prod(1);
        top(tc).V(i,1)=prod(2);
        top(tc).W(i,1)=prod(3);
        top(tc).P(i,1)=F4(i,j);

        top(tc).UNR(i,1)=F1(i,j);
        top(tc).VNR(i,1)=F2(i,j);
        top(tc).WNR(i,1)=F3(i,j);
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[F1(i,j); F2(i,j); F3(i,j)];
        bottom(bc).U(i,1)=prod(1);
        bottom(bc).V(i,1)=prod(2);
        bottom(bc).W(i,1)=prod(3);
        bottom(bc).P(i,1)=F4(i,j);

        bottom(bc).UNR(i,1)=F1(i,j);
        bottom(bc).VNR(i,1)=F2(i,j);
        bottom(bc).WNR(i,1)=F3(i,j);
    end
    
  end
  
% Pressure Gradient
  tc=0;
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[D7(i,j); D8(i,j); D9(i,j)];
        top(tc).dPdx(i,1)=prod(1);
        top(tc).dPdy(i,1)=prod(2);
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[D7(i,j); D8(i,j); D9(i,j)];
        bottom(bc).dPdx(i,1)=prod(1);
        bottom(bc).dPdy(i,1)=prod(2);
    end
    
  end

 % Laplace
  laplace_u = D45 + D46 + 0;
  laplace_v = D47 + D48 + 0;
  laplace_w = D49 + D50 + 0;
  tc=0;
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[laplace_u(i,j); laplace_v(i,j); laplace_w(i,j)];
        top(tc).LaplaceU(i,1)=prod(1);
        top(tc).LaplaceV(i,1)=prod(2);
        top(tc).LaplaceW(i,1)=prod(3);
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[laplace_u(i,j); laplace_v(i,j); laplace_w(i,j)];
        bottom(bc).LaplaceU(i,1)=prod(1);
        bottom(bc).LaplaceV(i,1)=prod(2);
        bottom(bc).LaplaceW(i,1)=prod(3);
    end
    
  end
 

% UU, VV, WW, UV, PP
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph)^2, sin(alph)^2,  2*cos(alph)*sin(alph), 0; ...
             sin(alph)^2, cos(alph)^2, -2*cos(alph)*sin(alph), 0; ...
             0 0 0 1; ...
            -cos(alph)*sin(alph), cos(alph)*sin(alph), (cos(alph)^2-sin(alph)^2), 0];
    for i=1:ln
        prod=Ralphap*[F5(i,j); F6(i,j); F9(i,j); F7(i,j)];
        top(tc).UU(i,1)=prod(1);
        top(tc).VV(i,1)=prod(2);
        top(tc).WW(i,1)=prod(3);
        top(tc).UV(i,1)=prod(4);
        top(tc).PP(i,1)=F8(i,j);

        top(tc).UUNR(i,1)=F5(i,j);
        top(tc).VVNR(i,1)=F6(i,j);
        top(tc).UVNR(i,1)=F9(i,j);
        top(tc).WWNR(i,1)=F7(i,j);
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph)^2, sin(alph)^2,  2*cos(alph)*sin(alph), 0; ...
             sin(alph)^2, cos(alph)^2, -2*cos(alph)*sin(alph), 0; ...
             0 0 0 1; ...
             -cos(alph)*sin(alph), cos(alph)*sin(alph), (cos(alph)^2-sin(alph)^2), 0];
    for i=1:ln
        prod=Ralphap*[F5(i,j); F6(i,j); F9(i,j); F7(i,j)];
        bottom(bc).UU(i,1)=prod(1);
        bottom(bc).VV(i,1)=prod(2);
        bottom(bc).WW(i,1)=prod(3);
        bottom(bc).UV(i,1)=prod(4);
        bottom(bc).PP(i,1)=F8(i,j);

        bottom(bc).UUNR(i,1)=F5(i,j);
        bottom(bc).VVNR(i,1)=F6(i,j);
        bottom(bc).UVNR(i,1)=F9(i,j);
        bottom(bc).WWNR(i,1)=F7(i,j);
    end
    
  end
  
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[F1(i,j); F2(i,j); F3(i,j)];
  %            bottom(j).U(i,1)=prod(1);
  %            bottom(j).V(i,1)=prod(2);
  %            bottom(j).W(i,1)=prod(3);
  %        end
  %end
  
  %Reynolds stress tensor. Tensors of Rank 2.
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                      (F9(i,j)-F1(i,j).*F2(i,j)) ...
                      (F11(i,j)-F1(i,j).*F3(i,j)); 
                      (F9(i,j)-F1(i,j).*F2(i,j)) ...
                      (F6(i,j)-F2(i,j).*F2(i,j)) ...
                      (F10(i,j)-F2(i,j).*F3(i,j)); 
                      (F11(i,j)-F1(i,j).*F3(i,j)) ...
                      (F10(i,j)-F2(i,j).*F3(i,j)) ...
                      (F7(i,j)-F3(i,j).*F3(i,j))]*Ralphap';
        top(tc).uu(i,1)=prod(1,1);
        top(tc).vv(i,1)=prod(2,2);
        top(tc).ww(i,1)=prod(3,3);
        top(tc).uv(i,1)=prod(1,2);
        top(tc).uw(i,1)=prod(1,3);
        top(tc).vw(i,1)=prod(2,3);

        prod=[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                      (F9(i,j)-F1(i,j).*F2(i,j)) ...
                      (F11(i,j)-F1(i,j).*F3(i,j)); 
                      (F9(i,j)-F1(i,j).*F2(i,j)) ...
                      (F6(i,j)-F2(i,j).*F2(i,j)) ...
                      (F10(i,j)-F2(i,j).*F3(i,j)); 
                      (F11(i,j)-F1(i,j).*F3(i,j)) ...
                      (F10(i,j)-F2(i,j).*F3(i,j)) ...
                      (F7(i,j)-F3(i,j).*F3(i,j))];
        top(tc).uuNR(i,1)=prod(1,1);
        top(tc).vvNR(i,1)=prod(2,2);
        top(tc).wwNR(i,1)=prod(3,3);
        top(tc).uvNR(i,1)=prod(1,2);
        top(tc).uwNR(i,1)=prod(1,3);
        top(tc).vwNR(i,1)=prod(2,3);


    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;    
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                       (F9(i,j)-F1(i,j).*F2(i,j)) ...
                       (F11(i,j)-F1(i,j).*F3(i,j)); 
                       (F9(i,j)-F1(i,j).*F2(i,j)) ...
                       (F6(i,j)-F2(i,j).*F2(i,j)) ...
                       (F10(i,j)-F2(i,j).*F3(i,j)); 
                       (F11(i,j)-F1(i,j).*F3(i,j)) ...
                       (F10(i,j)-F2(i,j).*F3(i,j)) ...
                       (F7(i,j)-F3(i,j).*F3(i,j))]*Ralphapp';
        bottom(bc).uu(i,1)=prod(1,1);
        bottom(bc).vv(i,1)=prod(2,2);
        bottom(bc).ww(i,1)=prod(3,3);
        bottom(bc).uv(i,1)=prod(1,2);
        bottom(bc).uw(i,1)=prod(1,3);
        bottom(bc).vw(i,1)=prod(2,3);

        prod=[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                       (F9(i,j)-F1(i,j).*F2(i,j)) ...
                       (F11(i,j)-F1(i,j).*F3(i,j)); 
                       (F9(i,j)-F1(i,j).*F2(i,j)) ...
                       (F6(i,j)-F2(i,j).*F2(i,j)) ...
                       (F10(i,j)-F2(i,j).*F3(i,j)); 
                       (F11(i,j)-F1(i,j).*F3(i,j)) ...
                       (F10(i,j)-F2(i,j).*F3(i,j)) ...
                       (F7(i,j)-F3(i,j).*F3(i,j))];
        bottom(bc).uuNR(i,1)=prod(1,1);
        bottom(bc).vvNR(i,1)=prod(2,2);
        bottom(bc).wwNR(i,1)=prod(3,3);
        bottom(bc).uvNR(i,1)=prod(1,2);
        bottom(bc).uwNR(i,1)=prod(1,3);
        bottom(bc).vwNR(i,1)=prod(2,3);

    end
    
  end

  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
  %                           (F9(i,j)-F1(i,j).*F2(i,j)) ...
  %                           (F11(i,j)-F1(i,j).*F3(i,j)); 
  %                           (F9(i,j)-F1(i,j).*F2(i,j)) ...
  %                           (F6(i,j)-F2(i,j).*F2(i,j)) ...
  %                           (F10(i,j)-F2(i,j).*F3(i,j)); 
  %                           (F11(i,j)-F1(i,j).*F3(i,j)) ...
  %                           (F10(i,j)-F2(i,j).*F3(i,j)) ...
  %                           (F7(i,j)-F3(i,j).*F3(i,j))]*Ralphapp';
  %            bottom(j).uu(i,1)=prod(1,1);
  %            bottom(j).vv(i,1)=prod(2,2);
  %            bottom(j).ww(i,1)=prod(3,3);
  %            bottom(j).uv(i,1)=prod(1,2);
  %            bottom(j).uw(i,1)=prod(1,3);
  %            bottom(j).vw(i,1)=prod(2,3);
  %        end
  %end
  
  %Mean, RMS, skewness and flatness of pressure
  P=F4;
  pp=F8-P.*P;
  ppp=F27-3*P.*pp-P.*P.*P;
  pppp=F38-4*P.*ppp-6*P.*P.*pp-P.*P.*P.*P;
  
  %Normalize pressure
  prms=sqrt(pp);
  pskew=ppp./(pp).^(3/2);
  pflat=pppp./(pp).^(2);
  
  %Velocity gradient tensor. Tensor of Rank 2.
  dUdx=D1;
  dVdx=D3;
  dWdx=D5;
  
  dUdy=D2;
  dVdy=D4;
  dWdy=D6;
  
  dUdz=zeros(size(dUdx));
  dVdz=zeros(size(dUdx));
  dWdz=zeros(size(dUdx));
  
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                      dVdx(i,j) dVdy(i,j) dVdz(i,j);
                      dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralphap';
                                        
        top(tc).dUdx(i,1)=prod(1,1);
        top(tc).dUdy(i,1)=prod(1,2);
        top(tc).dVdx(i,1)=prod(2,1);
        top(tc).dVdy(i,1)=prod(2,2);
        top(tc).dWdx(i,1)=prod(3,1);
        top(tc).dWdy(i,1)=prod(3,2);
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                       dVdx(i,j) dVdy(i,j) dVdz(i,j);
                       dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralphapp';
        
        
        bottom(bc).dUdx(i,1)=prod(1,1);
        bottom(bc).dUdy(i,1)=prod(1,2);
        bottom(bc).dVdx(i,1)=prod(2,1);
        bottom(bc).dVdy(i,1)=prod(2,2);
        bottom(bc).dWdx(i,1)=prod(3,1);
        bottom(bc).dWdy(i,1)=prod(3,2);
    end
    
  end


  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
  %                           dVdx(i,j) dVdy(i,j) dVdz(i,j);
  %                           dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralphapp';
  %
  %                      
  %            bottom(j).dUdx(i,1)=prod(1,1);
  %            bottom(j).dUdy(i,1)=prod(1,2);
  %            bottom(j).dVdx(i,1)=prod(2,1);
  %            bottom(j).dVdy(i,1)=prod(2,2);
  %            bottom(j).dWdx(i,1)=prod(3,1);
  %            bottom(j).dWdy(i,1)=prod(3,2);
  %        end
  %end
  
  %Production tensor. Tensor of Rank 2.
  %Definition of the production tensor assuming fully-developed flow, i.e.,
  %d()/dz=0, as a function of x and y.
  %Pij=-(<uiuk>*dUj/dxk+<ujuk>*dUi/dxk)
  %P11=-2*(<uu>*dU/dx+<uv>*dU/dy) 
  %P12=-(<uu>*dV/dx +<uv>*dV/dy+<uv>*dU/dx+<vv>*dU/dy)
  %P13=-(<uu>*dW/dx +<uv>*dW/dy+<uw>*dU/dx+<vw>*dU/dy)
  %P22=-2*(<uv>*dV/dx+<vv>*dV/dy) 
  %P23=-(<uv>*dW/dx+<vv>*dW/dy+<uw>*dV/dx+<vw>*dV/dy)
  %P33=-2*(<uw>*dW/dx+<vw>*dW/dy)
  
  %Reynolds stress tensor. Tensor of Rank 2.
  R2_tensor(1:3,1:3,1:ln,np)=0;
  for j=1:np
      for i=1:ln
          R2_tensor(:,:,i,j)=[ (F5(i,j)-F1(i,j).*F1(i,j)) (F9(i,j)-F1(i,j).*F2(i,j)) (F11(i,j)-F1(i,j).*F3(i,j));
                               (F9(i,j)-F1(i,j).*F2(i,j)) (F6(i,j)-F2(i,j).*F2(i,j)) (F10(i,j)-F2(i,j).*F3(i,j)); 
                              (F11(i,j)-F1(i,j).*F3(i,j)) (F10(i,j)-F2(i,j).*F3(i,j)) (F7(i,j)-F3(i,j).*F3(i,j))];
      end    
  end
  
  uu=squeeze(squeeze(R2_tensor(1,1,:,:)));
  vv=squeeze(squeeze(R2_tensor(2,2,:,:)));
  ww=squeeze(squeeze(R2_tensor(3,3,:,:)));
  uv=squeeze(squeeze(R2_tensor(1,2,:,:)));
  uw=squeeze(squeeze(R2_tensor(1,3,:,:)));
  vw=squeeze(squeeze(R2_tensor(2,3,:,:)));
  
  
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[-2*(uu(i,j).*dUdx(i,j)+uv(i,j).*dUdy(i,j)) ...
                       -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                       -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j));
                       -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                       -2*(uv(i,j).*dVdx(i,j)+vv(i,j).*dVdy(i,j)) ...
                       -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j));
                       -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j)) ...
                       -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j)) ...
                       -2*(uw(i,j).*dWdx(i,j)+vw(i,j).*dWdy(i,j))]*Ralphap';
                                        
        top(tc).Pxx(i,1)=prod(1,1);
        top(tc).Pyy(i,1)=prod(2,2);
        top(tc).Pzz(i,1)=prod(3,3);
        top(tc).Pxy(i,1)=prod(1,2);
        top(tc).Pxz(i,1)=prod(1,3);
        top(tc).Pyz(i,1)=prod(2,3);
  
    end
  end
  
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[-2*(uu(i,j).*dUdx(i,j)+uv(i,j).*dUdy(i,j)) ...
                       -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                       -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j));
                       -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                       -2*(uv(i,j).*dVdx(i,j)+vv(i,j).*dVdy(i,j)) ...
                       -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j));
                       -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j)) ...
                       -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j)) ...
                       -2*(uw(i,j).*dWdx(i,j)+vw(i,j).*dWdy(i,j))]*Ralphapp';
                    
        bottom(bc).Pxx(i,1)=prod(1,1);
        bottom(bc).Pyy(i,1)=prod(2,2);
        bottom(bc).Pzz(i,1)=prod(3,3);
        bottom(bc).Pxy(i,1)=prod(1,2);
        bottom(bc).Pxz(i,1)=prod(1,3);
        bottom(bc).Pyz(i,1)=prod(2,3);
  
    end
    
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[-2*(uu(i,j).*dUdx(i,j)+uv(i,j).*dUdy(i,j)) ...
  %                           -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
  %                           -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j));
  %                           -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
  %                           -2*(uv(i,j).*dVdx(i,j)+vv(i,j).*dVdy(i,j)) ...
  %                           -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j));
  %                           -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j)) ...
  %                           -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j)) ...
  %                           -2*(uw(i,j).*dWdx(i,j)+vw(i,j).*dWdy(i,j))]*Ralphapp';
  %
  %            bottom(j).Pxx(i,1)=prod(1,1);
  %            bottom(j).Pyy(i,1)=prod(2,2);
  %            bottom(j).Pzz(i,1)=prod(3,3);
  %            bottom(j).Pxy(i,1)=prod(1,2);
  %            bottom(j).Pxz(i,1)=prod(1,3);
  %            bottom(j).Pyz(i,1)=prod(2,3);
  %        end
  %end
  
  %Dissipation tensor. Tensor of Rank 2.
  %Definition of the dissipation tensor assuming fully-developed flow, i.e.,
  %d()/dz=0, as a function of x and y.
  %Dij=-2*nu*<dui/dxk*duj/dxk>
  
  %e11_tot=<((du/dx)^2+(du/dy)^2+(du/dz)^2)> % F42 
  %e22_tot=<((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)> % F43 
  %e33_tot=<((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)> % F44 
  %e12_tot=<(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F45  
  %e13_tot=<(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F46 
  %e23_tot=<(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F47 
  
  %e11=e11_tot-(dU/dx)^2-(dU/dy)^2 
  %e22=e22_tot-(dV/dx)^2-(dV/dy)^2 
  %e33=e33_tot-(dW/dx)^2-(dW/dy)^2 
  %e12=e12_tot-(dU/dx*dV/x)-(dU/dy*dV/dy)
  %e13=e13_tot-(dU/dx*dW/x)-(dU/dy*dW/dy)
  %e23=e23_tot-(dV/dx*dW/x)-(dV/dy*dW/dy)
  
  e11(1:ln,1:np)=0;
  e22(1:ln,1:np)=0;
  e33(1:ln,1:np)=0;
  e12(1:ln,1:np)=0;
  e13(1:ln,1:np)=0;
  e23(1:ln,1:np)=0;
  
  for j=1:np
      for i=1:ln
          dUidxj(:,:,i,j)=[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                           dVdx(i,j) dVdy(i,j) dVdz(i,j);
                           dWdx(i,j) dWdy(i,j) dWdz(i,j)];
      end
  end
  
  for j=1:np
      e11(:,j)=F42(:,j)- ...
          (squeeze(dUidxj(1,1,:,j))).^2-(squeeze(dUidxj(1,2,:,j))).^2; 
      e22(:,j)=F43(:,j)- ...
          (squeeze(dUidxj(2,1,:,j))).^2-(squeeze(dUidxj(2,2,:,j))).^2; 
      e33(:,j)=F44(:,j)- ...
          (squeeze(dUidxj(3,1,:,j))).^2-(squeeze(dUidxj(3,2,:,j))).^2; 
      e12(:,j)=F45(:,j)- ...
          squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(2,1,:,j))- ...
          squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(2,2,:,j));
      e13(:,j)=F46(:,j)- ...
          squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(3,1,:,j))- ...
          squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(3,2,:,j));
      e23(:,j)=F47(:,j)- ...
          squeeze(dUidxj(2,1,:,j)).*squeeze(dUidxj(3,1,:,j))- ...
          squeeze(dUidxj(2,2,:,j)).*squeeze(dUidxj(3,2,:,j));      
  end
  
  for tc=1:topcount 
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[e11(i,j) e12(i,j) e13(i,j);
                      e12(i,j) e22(i,j) e23(i,j);
                      e13(i,j) e23(i,j) e33(i,j)]*Ralphap';
                                        
        top(tc).Dxx(i,1)=-2*nu*prod(1,1);
        top(tc).Dyy(i,1)=-2*nu*prod(2,2);
        top(tc).Dzz(i,1)=-2*nu*prod(3,3);
        top(tc).Dxy(i,1)=-2*nu*prod(1,2);
        top(tc).Dxz(i,1)=-2*nu*prod(1,3);
        top(tc).Dyz(i,1)=-2*nu*prod(2,3);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[e11(i,j) e12(i,j) e13(i,j);
                       e12(i,j) e22(i,j) e23(i,j);
                       e13(i,j) e23(i,j) e33(i,j)]*Ralphapp';
                    
        bottom(bc).Dxx(i,1)=-2*nu*prod(1,1);
        bottom(bc).Dyy(i,1)=-2*nu*prod(2,2);
        bottom(bc).Dzz(i,1)=-2*nu*prod(3,3);
        bottom(bc).Dxy(i,1)=-2*nu*prod(1,2);
        bottom(bc).Dxz(i,1)=-2*nu*prod(1,3);
        bottom(bc).Dyz(i,1)=-2*nu*prod(2,3);
  
    end
    
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[e11(i,j) e12(i,j) e13(i,j);
  %                           e12(i,j) e22(i,j) e23(i,j);
  %                           e13(i,j) e23(i,j) e33(i,j)]*Ralphapp';
  %
  %            bottom(j).Dxx(i,1)=-2*nu*prod(1,1);
  %            bottom(j).Dyy(i,1)=-2*nu*prod(2,2);
  %            bottom(j).Dzz(i,1)=-2*nu*prod(3,3);
  %            bottom(j).Dxy(i,1)=-2*nu*prod(1,2);
  %            bottom(j).Dxz(i,1)=-2*nu*prod(1,3);
  %            bottom(j).Dyz(i,1)=-2*nu*prod(2,3);
  %        end
  %end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Additional dissipation from the RT filter
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Dissipation from RT filter. Tensors of Rank 1.
  if (nstat==66)
    for tc=1:topcount 
      j=tindicies(tc);alph=top(tc).alpha;
      Ralphap=[cos(alph) sin(alph) 0; ...
          -sin(alph) cos(alph) 0; ...
          0 0 1];
      for i=1:ln
          prod=Ralphap*[F64(i,j); F65(i,j); F66(i,j)];
          top(tc).DRTx(i,1)=prod(1);
          top(tc).DRTy(i,1)=prod(2);
          top(tc).DRTz(i,1)=prod(3);
      end
    end
    
    for bc=1:botcount
      j=bindicies(bc);alph=bottom(bc).alpha;
      Ralphapp=[cos(alph) sin(alph) 0; ...
          sin(alph) -cos(alph) 0; ...
          0 0 1];
      for i=1:ln
          prod=Ralphapp*[F64(i,j); F65(i,j); F66(i,j)];
          bottom(bc).DRTx(i,1)=prod(1);
          bottom(bc).DRTy(i,1)=prod(2);
          bottom(bc).DRTz(i,1)=prod(3);
      end
    end
  end       % nstat==66  
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[F64(i,j); F65(i,j); F66(i,j)];
  %            bottom(j).DRTx(i,1)=prod(1);
  %            bottom(j).DRTy(i,1)=prod(2);
  %            bottom(j).DRTz(i,1)=prod(3);
  %        end
  %end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %Mean convection tensor. Tensor of Rank 2.
  %Definition of the mean convection tensor assuming 
  %fully-developed flow, i.e., d()/dz=0, as a function of x and y.
  %Cij=Uk*d<uiuj>/dxk
  %Note that under this definition: Production + Dissipation - Convection ...
  %C11=U*d(uu)/dx+V*d(uu)/dy
  %C22=U*d(vv)/dx+V*d(vv)/dy
  %C33=U*d(ww)/dx+V*d(ww)/dy
  %C12=U*d(uv)/dx+V*d(uv)/dy
  %C13=U*d(uw)/dx+V*d(uw)/dy
  %C23=U*d(vw)/dx+V*d(vw)/dy
  
  U=F1;
  V=F2;
  W=F3;
  
  duudx=D9-2*U.*dUdx;
  dvvdx=D11-2*V.*dVdx;
  dwwdx=D13-2*W.*dWdx;
  duvdx=D17-U.*dVdx-V.*dUdx;
  duwdx=D21-U.*dWdx-W.*dUdx;
  dvwdx=D19-V.*dWdx-W.*dVdx;
  
  duudy=D10-2*U.*dUdy;
  dvvdy=D12-2*V.*dVdy;
  dwwdy=D14-2*W.*dWdy;
  duvdy=D18-U.*dVdy-V.*dUdy;
  duwdy=D22-U.*dWdy-W.*dUdy;
  dvwdy=D20-V.*dWdy-W.*dVdy;
  
  for tc=1:topcount 
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[duudx(i,j) duudy(i,j) 0; ...
                      dvvdx(i,j) dvvdy(i,j) 0; ...
                      dwwdx(i,j) dwwdy(i,j) 0; ]*Ralphap';
        top(tc).duudx(i,1)=prod(1,1);
        top(tc).duudy(i,1)=prod(1,2);
        top(tc).dvvdx(i,1)=prod(2,1);
        top(tc).dvvdy(i,1)=prod(2,2);
        top(tc).dwwdx(i,1)=prod(3,1);
        top(tc).dwwdy(i,1)=prod(3,2);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[duudx(i,j) duudy(i,j) 0; ...
                      dvvdx(i,j) dvvdy(i,j) 0; ...
                      dwwdx(i,j) dwwdy(i,j) 0; ]*Ralphap';
        bottom(bc).duudx(i,1)=prod(1,1);
        bottom(bc).duudy(i,1)=prod(1,2);
        bottom(bc).dvvdx(i,1)=prod(2,1);
        bottom(bc).dvvdy(i,1)=prod(2,2);
        bottom(bc).dwwdx(i,1)=prod(3,1);
        bottom(bc).dwwdy(i,1)=prod(3,2);
  
    end
  end
  
  
  duudx_les=duudx;
  dvvdx_les=dvvdx;
  dwwdx_les=dwwdx;
  duvdx_les=duvdx;
  duwdx_les=duwdx;
  dvwdx_les=dvwdx;
  
  duudy_les=duudy;
  dvvdy_les=dvvdy;
  dwwdy_les=dwwdy;
  duvdy_les=duvdy;
  duwdy_les=duwdy;
  dvwdy_les=dvwdy;
  
  U_les=U;
  V_les=V;
  
  duudz=zeros(size(duudx));
  dvvdz=zeros(size(duudx));
  dwwdz=zeros(size(duudx));
  duvdz=zeros(size(duudx));
  duwdz=zeros(size(duudx));
  dvwdz=zeros(size(duudx));
  
  for tc=1:topcount 
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[U(i,j)*duudx(i,j)+V(i,j)*duudy(i,j) ...
                      U(i,j)*duvdx(i,j)+V(i,j)*duvdy(i,j) ...
                      U(i,j)*duwdx(i,j)+V(i,j)*duwdy(i,j);
                      U(i,j)*duvdx(i,j)+V(i,j)*duvdy(i,j) ...
                      U(i,j)*dvvdx(i,j)+V(i,j)*dvvdy(i,j) ...
                      U(i,j)*dvwdx(i,j)+V(i,j)*dvwdy(i,j);
                      U(i,j)*duwdx(i,j)+V(i,j)*duwdy(i,j) ...
                      U(i,j)*dvwdx(i,j)+V(i,j)*dvwdy(i,j) ...
                      U(i,j)*dwwdx(i,j)+V(i,j)*dwwdy(i,j)]*Ralphap';
                                        
        top(tc).Cxx(i,1)=prod(1,1);
        top(tc).Cyy(i,1)=prod(2,2);
        top(tc).Czz(i,1)=prod(3,3);
        top(tc).Cxy(i,1)=prod(1,2);
        top(tc).Cxz(i,1)=prod(1,3);
        top(tc).Cyz(i,1)=prod(2,3);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[U(i,j)*duudx(i,j)+V(i,j)*duudy(i,j) ...
                       U(i,j)*duvdx(i,j)+V(i,j)*duvdy(i,j) ...
                       U(i,j)*duwdx(i,j)+V(i,j)*duwdy(i,j);
                       U(i,j)*duvdx(i,j)+V(i,j)*duvdy(i,j) ...
                       U(i,j)*dvvdx(i,j)+V(i,j)*dvvdy(i,j) ...
                       U(i,j)*dvwdx(i,j)+V(i,j)*dvwdy(i,j);
                       U(i,j)*duwdx(i,j)+V(i,j)*duwdy(i,j) ...
                       U(i,j)*dvwdx(i,j)+V(i,j)*dvwdy(i,j) ...
                       U(i,j)*dwwdx(i,j)+V(i,j)*dwwdy(i,j)]*Ralphapp';
                    
        bottom(bc).Cxx(i,1)=prod(1,1);
        bottom(bc).Cyy(i,1)=prod(2,2);
        bottom(bc).Czz(i,1)=prod(3,3);
        bottom(bc).Cxy(i,1)=prod(1,2);
        bottom(bc).Cxz(i,1)=prod(1,3);
        bottom(bc).Cyz(i,1)=prod(2,3);
  
    end
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[U(i,j)*duudx(i,j)+V(i,j)*duudy(i,j) ...
  %                           U(i,j)*duvdx(i,j)+V(i,j)*duvdy(i,j) ...
  %                           U(i,j)*duwdx(i,j)+V(i,j)*duwdy(i,j);
  %                           U(i,j)*duvdx(i,j)+V(i,j)*duvdy(i,j) ...
  %                           U(i,j)*dvvdx(i,j)+V(i,j)*dvvdy(i,j) ...
  %                           U(i,j)*dvwdx(i,j)+V(i,j)*dvwdy(i,j);
  %                           U(i,j)*duwdx(i,j)+V(i,j)*duwdy(i,j) ...
  %                           U(i,j)*dvwdx(i,j)+V(i,j)*dvwdy(i,j) ...
  %                           U(i,j)*dwwdx(i,j)+V(i,j)*dwwdy(i,j)]*Ralphapp';
  %
  %            bottom(j).Cxx(i,1)=prod(1,1);
  %            bottom(j).Cyy(i,1)=prod(2,2);
  %            bottom(j).Czz(i,1)=prod(3,3);
  %            bottom(j).Cxy(i,1)=prod(1,2);
  %            bottom(j).Cxz(i,1)=prod(1,3);
  %            bottom(j).Cyz(i,1)=prod(2,3);
  %        end
  %end
  
  %Turbulent transport tensor. Tensor of Rank 2.
  %Definition of the turbulent transport tensor assuming fully-developed 
  %flow, i.e., d()/dz=0, as a function of x and y.
  %Tij=-d<uiujuk>/dxk
  %T11=-(d<uuuk>)/dxk, K=1,3=-[d<uuu>/dx+d<uuv>/dy+d<uuw>/dz]
  %T22=-(d<vvuk>)/dxk, K=1,3=-[d<vvu>/dx+d<vvv>/dy+d<vvw>/dz]
  %T33=-(d<wwuk>)/dxk, K=1,3=-[d<wwu>/dx+d<wwv>/dy+d<www>/dz]
  %T12=-(d<uvuk>)/dxk, K=1,3=-[d<uvu>/dx+d<uvv>/dy+d<uvw>/dz]
  %T13=-(d<uwuk>)/dxk, K=1,3=-[d<uwu>/dx+d<uwv>/dy+d<uww>/dz]
  %T23=-(d<vwuk>)/dxk, K=1,3=-[d<vwu>/dx+d<vwv>/dy+d<vww>/dz]
  
  duuudx=D23-3*U.*U.*dUdx-3*(U.*duudx+uu.*dUdx);
  dvvudx=D35-2*(V.*duvdx+uv.*dVdx)-(U.*dvvdx+vv.*dUdx) ...
      -(V.*V.*dUdx+2*U.*V.*dVdx);
  dwwudx=D39-2*(W.*duwdx+uw.*dWdx)-(U.*dwwdx+ww.*dUdx) ... 
      -(W.*W.*dUdx+2*U.*W.*dWdx);
  duvudx=D31-2*(U.*duvdx+uv.*dUdx)-(V.*duudx+uu.*dVdx) ...
      -(U.*U.*dVdx+2*U.*V.*dUdx);
  duwudx=D33-2*(U.*duwdx+uw.*dUdx)-(W.*duudx+uu.*dWdx) ...
      -(U.*U.*dWdx+2*U.*W.*dUdx);
  dvwudx=D43-(U.*dvwdx+vw.*dUdx)-(V.*duwdx+uw.*dVdx)-(W.*duvdx+uv.*dWdx) ...
      -(U.*V.*dWdx+U.*W.*dVdx+V.*W.*dUdx);
        
  duuvdy=D32-2*(U.*duvdy+uv.*dUdy)-(V.*duudy+uu.*dVdy) ...
      -(U.*U.*dVdy+2*U.*V.*dUdy);
  dvvvdy=D26-3*(V.*dvvdy+vv.*dVdy)-3*V.*V.*dVdy;
  dwwvdy=D42-2*(W.*dvwdy+vw.*dWdy)-(V.*dwwdy+ww.*dVdy) ...
      -(W.*W.*dVdy+2*V.*W.*dWdy);
  duvvdy=D36-2*(V.*duvdy+uv.*dVdy)-(U.*dvvdy+vv.*dUdy) ...
      -(V.*V.*dUdy+2*U.*V.*dVdy);
  duwvdy=D44-(U.*dvwdy+vw.*dUdy)-(V.*duwdy+uw.*dVdy)-(W.*duvdy+uv.*dWdy) ...
      -(U.*V.*dWdy+U.*W.*dVdy+V.*W.*dUdy);
  dvwvdy=D38-2*(V.*dvwdy+vw.*dVdy)-(W.*dvvdy+vv.*dWdy) ...
      -(V.*V.*dWdy+2*V.*W.*dVdy);
  
  duuwdz=zeros(size(duuudx));
  dvvwdz=zeros(size(duuudx));
  dwwwdz=zeros(size(duuudx));
  duvwdz=zeros(size(duuudx));
  duwwdz=zeros(size(duuudx));
  dvwwdz=zeros(size(duuudx));
  
  duuudx_les=duuudx;
  dvvvdy_les=dvvvdy;
  
  for tc=1:topcount 
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[-(duuudx(i,j)+duuvdy(i,j)+duuwdz(i,j)) ...
                      -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                      -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j));
                      -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                      -(dvvudx(i,j)+dvvvdy(i,j)+dvvwdz(i,j)) ...
                      -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j));
                      -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j)) ...
                      -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j)) ...
                      -(dwwudx(i,j)+dwwvdy(i,j)+dwwwdz(i,j))]*Ralphap';
                                        
        top(tc).Txx(i,1)=prod(1,1);
        top(tc).Tyy(i,1)=prod(2,2);
        top(tc).Tzz(i,1)=prod(3,3);
        top(tc).Txy(i,1)=prod(1,2);
        top(tc).Txz(i,1)=prod(1,3);
        top(tc).Tyz(i,1)=prod(2,3);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[-(duuudx(i,j)+duuvdy(i,j)+duuwdz(i,j)) ...
                       -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                       -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j));
                       -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                       -(dvvudx(i,j)+dvvvdy(i,j)+dvvwdz(i,j)) ...
                       -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j));
                       -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j)) ...
                       -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j)) ...
                       -(dwwudx(i,j)+dwwvdy(i,j)+dwwwdz(i,j))]*Ralphapp';
                    
        bottom(bc).Txx(i,1)=prod(1,1);
        bottom(bc).Tyy(i,1)=prod(2,2);
        bottom(bc).Tzz(i,1)=prod(3,3);
        bottom(bc).Txy(i,1)=prod(1,2);
        bottom(bc).Txz(i,1)=prod(1,3);
        bottom(bc).Tyz(i,1)=prod(2,3);
    end
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[-(duuudx(i,j)+duuvdy(i,j)+duuwdz(i,j)) ...
  %                           -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
  %                           -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j));
  %                           -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
  %                           -(dvvudx(i,j)+dvvvdy(i,j)+dvvwdz(i,j)) ...
  %                           -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j));
  %                           -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j)) ...
  %                           -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j)) ...
  %                           -(dwwudx(i,j)+dwwvdy(i,j)+dwwwdz(i,j))]*Ralphapp';
  %
  %            bottom(j).Txx(i,1)=prod(1,1);
  %            bottom(j).Tyy(i,1)=prod(2,2);
  %            bottom(j).Tzz(i,1)=prod(3,3);
  %            bottom(j).Txy(i,1)=prod(1,2);
  %            bottom(j).Txz(i,1)=prod(1,3);
  %            bottom(j).Tyz(i,1)=prod(2,3);
  %        end
  %end
  
  %Viscous diffusion tensor. Tensor of Rank 2.
  %Definition of the viscous diffusion tensor assuming fully-developed 
  %flow, i.e., d()/dz=0, as a function of x and y.
  %VDij=nu*d2(uiuj)/dxk2
  %VD11=nu*(d2(uu)/dx2+d2(uu)/dy2+d2(uu)/dz2)
  %VD22=nu*(d2(vv)/dx2+d2(vv)/dy2+d2(vv)/dz2)
  %VD33=nu*(d2(ww)/dx2+d2(ww)/dy2+d2(ww)/dz2)
  %VD12=nu*(d2(uv)/dx2+d2(uv)/dy2+d2(uv)/dz2)
  %VD13=nu*(d2(uw)/dx2+d2(uw)/dy2+d2(uw)/dz2)
  %VD23=nu*(d2(vw)/dx2+d2(vw)/dy2+d2(vw)/dz2)
  
  d2uudx2=D51-2*(U.*D45+dUdx.*dUdx);
  d2vvdx2=D53-2*(V.*D47+dVdx.*dVdx);
  d2wwdx2=D55-2*(W.*D49+dWdx.*dWdx);
  d2uvdx2=D57-(V.*D45+U.*D47+2*dUdx.*dVdx);
  d2uwdx2=D59-(U.*D49+W.*D45+2*dUdx.*dWdx);
  d2vwdx2=D61-(V.*D49+W.*D47+2*dVdx.*dWdx);
  
  d2uudy2=D52-2*(U.*D46+dUdy.*dUdy);
  d2vvdy2=D54-2*(V.*D48+dVdy.*dVdy);
  d2wwdy2=D56-2*(W.*D50+dWdy.*dWdy);
  d2uvdy2=D58-(V.*D46+U.*D48+2*dUdy.*dVdy);
  d2uwdy2=D60-(U.*D50+W.*D46+2*dUdy.*dWdy);
  d2vwdy2=D62-(V.*D50+W.*D48+2*dVdy.*dWdy);
  
  d2uudz2=zeros(size(d2uudx2));
  d2vvdz2=zeros(size(d2uudx2));
  d2wwdz2=zeros(size(d2uudx2));
  d2uvdz2=zeros(size(d2uudx2));
  d2uwdz2=zeros(size(d2uudx2));
  d2vwdz2=zeros(size(d2uudx2));
  
  for tc=1:topcount 
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[nu*(d2uudx2(i,j)+d2uudy2(i,j)+d2uudz2(i,j)) ...
                      nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                      nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j));
                      nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                      nu*(d2vvdx2(i,j)+d2vvdy2(i,j)+d2vvdz2(i,j)) ...
                      nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j));
                      nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j)) ...
                      nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j)) ...
                      nu*(d2wwdx2(i,j)+d2wwdy2(i,j)+d2wwdz2(i,j))]*Ralphap';
                                        
        top(tc).VDxx(i,1)=prod(1,1);
        top(tc).VDyy(i,1)=prod(2,2);
        top(tc).VDzz(i,1)=prod(3,3);
        top(tc).VDxy(i,1)=prod(1,2);
        top(tc).VDxz(i,1)=prod(1,3);
        top(tc).VDyz(i,1)=prod(2,3);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[nu*(d2uudx2(i,j)+d2uudy2(i,j)+d2uudz2(i,j)) ...
                       nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                       nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j));
                       nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                       nu*(d2vvdx2(i,j)+d2vvdy2(i,j)+d2vvdz2(i,j)) ...
                       nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j));
                       nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j)) ...
                       nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j)) ...
                       nu*(d2wwdx2(i,j)+d2wwdy2(i,j)+d2wwdz2(i,j))]*Ralphapp';
                    
        bottom(bc).VDxx(i,1)=prod(1,1);
        bottom(bc).VDyy(i,1)=prod(2,2);
        bottom(bc).VDzz(i,1)=prod(3,3);
        bottom(bc).VDxy(i,1)=prod(1,2);
        bottom(bc).VDxz(i,1)=prod(1,3);
        bottom(bc).VDyz(i,1)=prod(2,3);
    end
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[nu*(d2uudx2(i,j)+d2uudy2(i,j)+d2uudz2(i,j)) ...
  %                           nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
  %                           nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j));
  %                           nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
  %                           nu*(d2vvdx2(i,j)+d2vvdy2(i,j)+d2vvdz2(i,j)) ...
  %                           nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j));
  %                           nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j)) ...
  %                           nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j)) ...
  %                           nu*(d2wwdx2(i,j)+d2wwdy2(i,j)+d2wwdz2(i,j))]*Ralphapp';
  %
  %            bottom(j).VDxx(i,1)=prod(1,1);
  %            bottom(j).VDyy(i,1)=prod(2,2);
  %            bottom(j).VDzz(i,1)=prod(3,3);
  %            bottom(j).VDxy(i,1)=prod(1,2);
  %            bottom(j).VDxz(i,1)=prod(1,3);
  %            bottom(j).VDyz(i,1)=prod(2,3);
  %        end
  %end
  
  %Velocity-pressure-gradient tensor. Tensor of Rank 2.
  %Definition of the velocity-pressure-gradient tensor assuming 
  %fully-developed flow, i.e., d()/dz=0, as a function of x and y.
  %Piij=-1/rho*(<ui*dp/dxj>+<uj*dp/dxi>)
  %Pi11=-2/rho*<u*dp/dx>
  %Pi22=-2/rho*<v*dp/dy>
  %Pi33=-2/rho*<w*dp/dz>
  %Pi12=-1/rho*(<u*dp/dy>+<v*dp/dx>)
  %Pi13=-1/rho*(<u*dp/dz>+<w*dp/dx>)
  %Pi23=-1/rho*(<v*dp/dz>+<w*dp/dy>)
  
  %Now, since we don't compute <ui*dp/dxj>, we use the chain rule to express
  %these terms as a function of velocity gradients.
  %<ui*dp/dxj>=d(<p*ui>)/dxj-<p*dui/dxj>
  
  %We define the pressure transport and pressure strain tensors.
  
  %Pressure transport tensor. Tensor of Rank 2.
  %Definition of the pressure transport tensor assuming 
  %fully-developed flow, i.e., d()/dz=0, as a function of x and y.
  %PTij=-1/rho*(<d(p*ui)/dxj>+<d(p*uj)/dxi>)
  %PT11=-2/rho*<d(p*u)/dx>
  %PT22=-2/rho*<d(p*v)/dy>
  %PT33=-2/rho*<d(p*w)/dz>
  %PT12=-1/rho*(<d(p*u)/dy>+<d(p*v)/dx>)
  %PT13=-1/rho*(<d(p*u)/dz>+<d(p*w)/dx>)
  %PT23=-1/rho*(<d(p*v)/dz>+<d(p*w)/dy>)
  
  dpudx=D63-P.*dUdx-U.*D7;
  dpvdx=D65-P.*dVdx-V.*D7;
  dpwdx=D67-P.*dWdx-W.*D7;
  
  dpudy=D64-P.*dUdy-U.*D8;
  dpvdy=D66-P.*dVdy-V.*D8;
  dpwdy=D68-P.*dWdy-W.*D8;
  
  dpudz=zeros(size(dpudx));
  dpvdz=zeros(size(dpudx));
  dpwdz=zeros(size(dpudx));
  
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[(dpudx(i,j)) ... 
                      (dpudy(i,j)+dpvdx(i,j)) ...
                      (dpudz(i,j)+dpwdx(i,j));
                      (dpudy(i,j)+dpvdx(i,j)) ...
                      (dpvdy(i,j)) ...
                      (dpvdz(i,j)+dpwdy(i,j));
                      (dpudz(i,j)+dpwdx(i,j)) ...
                      (dpvdz(i,j)+dpwdy(i,j)) ...
                      (dpwdz(i,j))]*Ralphap';
                                        
        top(tc).PTxx(i,1)=prod(1,1);
        top(tc).PTyy(i,1)=prod(2,2);
        top(tc).PTzz(i,1)=prod(3,3);
        top(tc).PTxy(i,1)=prod(1,2);
        top(tc).PTxz(i,1)=prod(1,3);
        top(tc).PTyz(i,1)=prod(2,3);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;  
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[(dpudx(i,j)) ... 
                       (dpudy(i,j)+dpvdx(i,j)) ...
                       (dpudz(i,j)+dpwdx(i,j));
                       (dpudy(i,j)+dpvdx(i,j)) ...
                       (dpvdy(i,j)) ...
                       (dpvdz(i,j)+dpwdy(i,j));
                       (dpudz(i,j)+dpwdx(i,j)) ...
                       (dpvdz(i,j)+dpwdy(i,j)) ...
                       (dpwdz(i,j))]*Ralphapp';
                    
        bottom(bc).PTxx(i,1)=prod(1,1);
        bottom(bc).PTyy(i,1)=prod(2,2);
        bottom(bc).PTzz(i,1)=prod(3,3);
        bottom(bc).PTxy(i,1)=prod(1,2);
        bottom(bc).PTxz(i,1)=prod(1,3);
        bottom(bc).PTyz(i,1)=prod(2,3);
  
    end
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[(dpudx(i,j)) ... 
  %                           (dpudy(i,j)+dpvdx(i,j)) ...
  %                           (dpudz(i,j)+dpwdx(i,j));
  %                           (dpudy(i,j)+dpvdx(i,j)) ...
  %                           (dpvdy(i,j)) ...
  %                           (dpvdz(i,j)+dpwdy(i,j));
  %                           (dpudz(i,j)+dpwdx(i,j)) ...
  %                           (dpvdz(i,j)+dpwdy(i,j)) ...
  %                           (dpwdz(i,j))]*Ralphapp';
  %
  %            bottom(j).PTxx(i,1)=prod(1,1);
  %            bottom(j).PTyy(i,1)=prod(2,2);
  %            bottom(j).PTzz(i,1)=prod(3,3);
  %            bottom(j).PTxy(i,1)=prod(1,2);
  %            bottom(j).PTxz(i,1)=prod(1,3);
  %            bottom(j).PTyz(i,1)=prod(2,3);
  %        end
  %end
  
  %Pressure strain tensor. Tensor of Rank 2.
  %Definition of the pressure strain tensor assuming 
  %fully-developed flow, i.e., d()/dz=0, as a function of x and y.
  %PSij=-1/rho*(<p*dui/dxj>+<p*duj/dxi>)
  %PS11=-2/rho*<p*du/dx>
  %PS22=-2/rho*<p*dv/dy>
  %PS33=-2/rho*<p*dw/dz>
  %PS12=-1/rho*(<p*du/dy>+<p*dv/dx>)
  %PS13=-1/rho*(<p*du/dz>+<p*dw/dx>)
  %PS23=-1/rho*(<p*dv/dz>+<p*dw/dy>)
  
  %<pdudx>       % F15
  %<pdudy>       % F16
  %<pdudz>       % F17
  
  %<pdvdx>       % F18
  %<pdvdy>       % F19
  %<pdvdz>       % F20
  
  %<pdwdx>       % F21
  %<pdwdy>       % F22
  %<pdwdz>       % F23
  
  pdudx=F15-P.*dUdx;
  pdudy=F16-P.*dUdy;
  pdudz=F17-P.*dUdz;
  
  pdvdx=F18-P.*dVdx;
  pdvdy=F19-P.*dVdy;
  pdvdz=F20-P.*dVdz;
  
  pdwdx=F21-P.*dWdx;
  pdwdy=F22-P.*dWdy;
  pdwdz=F23-P.*dWdz;
  
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphap*[(pdudx(i,j)) ... 
                      (pdudy(i,j)+pdvdx(i,j)) ...
                      (pdudz(i,j)+pdwdx(i,j));
                      (pdudy(i,j)+pdvdx(i,j)) ...
                      (pdvdy(i,j)) ...
                      (pdvdz(i,j)+pdwdy(i,j));
                      (pdudz(i,j)+pdwdx(i,j)) ...
                      (pdvdz(i,j)+pdwdy(i,j)) ...
                      (pdwdz(i,j))]*Ralphap';
                                        
        top(tc).PSxx(i,1)=prod(1,1);
        top(tc).PSyy(i,1)=prod(2,2);
        top(tc).PSzz(i,1)=prod(3,3);
        top(tc).PSxy(i,1)=prod(1,2);
        top(tc).PSxz(i,1)=prod(1,3);
        top(tc).PSyz(i,1)=prod(2,3);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        prod=Ralphapp*[(pdudx(i,j)) ... 
                       (pdudy(i,j)+pdvdx(i,j)) ...
                       (pdudz(i,j)+pdwdx(i,j));
                       (pdudy(i,j)+pdvdx(i,j)) ...
                       (pdvdy(i,j)) ...
                       (pdvdz(i,j)+pdwdy(i,j));
                       (pdudz(i,j)+pdwdx(i,j)) ...
                       (pdvdz(i,j)+pdwdy(i,j)) ...
                       (pdwdz(i,j))]*Ralphapp';
                    
        bottom(bc).PSxx(i,1)=prod(1,1);
        bottom(bc).PSyy(i,1)=prod(2,2);
        bottom(bc).PSzz(i,1)=prod(3,3);
        bottom(bc).PSxy(i,1)=prod(1,2);
        bottom(bc).PSxz(i,1)=prod(1,3);
        bottom(bc).PSyz(i,1)=prod(2,3);
  
    end
    
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            prod=Ralphapp*[(pdudx(i,j)) ... 
  %                           (pdudy(i,j)+pdvdx(i,j)) ...
  %                           (pdudz(i,j)+pdwdx(i,j));
  %                           (pdudy(i,j)+pdvdx(i,j)) ...
  %                           (pdvdy(i,j)) ...
  %                           (pdvdz(i,j)+pdwdy(i,j));
  %                           (pdudz(i,j)+pdwdx(i,j)) ...
  %                           (pdvdz(i,j)+pdwdy(i,j)) ...
  %                           (pdwdz(i,j))]*Ralphapp';
  %
  %            bottom(j).PSxx(i,1)=prod(1,1);
  %            bottom(j).PSyy(i,1)=prod(2,2);
  %            bottom(j).PSzz(i,1)=prod(3,3);
  %            bottom(j).PSxy(i,1)=prod(1,2);
  %            bottom(j).PSxz(i,1)=prod(1,3);
  %            bottom(j).PSyz(i,1)=prod(2,3);
  %        end
  %end
  
  %Construct velocity-pressure-gradient-tensor
  for i=1:length(top)
      top(i).Pixx=-2/rho*(top(i).PTxx-top(i).PSxx);
      top(i).Piyy=-2/rho*(top(i).PTyy-top(i).PSyy);
      top(i).Pizz=-2/rho*(top(i).PTzz-top(i).PSzz);
      top(i).Pixy=-1/rho*(top(i).PTxy-top(i).PSxy);
      top(i).Pixz=-1/rho*(top(i).PTxz-top(i).PSxz);
      top(i).Piyz=-1/rho*(top(i).PTyz-top(i).PSyz);
  end
  
  for i=1:length(bottom)    
      bottom(i).Pixx=-2/rho*(bottom(i).PTxx-bottom(i).PSxx);
      bottom(i).Piyy=-2/rho*(bottom(i).PTyy-bottom(i).PSyy);
      bottom(i).Pizz=-2/rho*(bottom(i).PTzz-bottom(i).PSzz);
      bottom(i).Pixy=-1/rho*(bottom(i).PTxy-bottom(i).PSxy);
      bottom(i).Pixz=-1/rho*(bottom(i).PTxz-bottom(i).PSxz);
      bottom(i).Piyz=-1/rho*(bottom(i).PTyz-bottom(i).PSyz);
  end
  
  %Budget for each component of the Reynolds stress tensor 
  %Without mean convection
  for i=1:length(top)
      top(i).Sxx=top(i).Pxx+top(i).Dxx+top(i).Txx+top(i).VDxx+top(i).Pixx; 
      top(i).Syy=top(i).Pyy+top(i).Dyy+top(i).Tyy+top(i).VDyy+top(i).Piyy;
      top(i).Szz=top(i).Pzz+top(i).Dzz+top(i).Tzz+top(i).VDzz+top(i).Pizz;
      top(i).Sxy=top(i).Pxy+top(i).Dxy+top(i).Txy+top(i).VDxy+top(i).Pixy; 
      top(i).Sxz=top(i).Pxz+top(i).Dxz+top(i).Txz+top(i).VDxz+top(i).Pixz;
      top(i).Syz=top(i).Pyz+top(i).Dyz+top(i).Tyz+top(i).VDyz+top(i).Piyz;
  end
  
  for i=1:length(bottom)    
      bottom(i).Sxx=bottom(i).Pxx+bottom(i).Dxx+bottom(i).Txx+bottom(i).VDxx+bottom(i).Pixx; 
      bottom(i).Syy=bottom(i).Pyy+bottom(i).Dyy+bottom(i).Tyy+bottom(i).VDyy+bottom(i).Piyy;
      bottom(i).Szz=bottom(i).Pzz+bottom(i).Dzz+bottom(i).Tzz+bottom(i).VDzz+bottom(i).Pizz;
      bottom(i).Sxy=bottom(i).Pxy+bottom(i).Dxy+bottom(i).Txy+bottom(i).VDxy+bottom(i).Pixy; 
      bottom(i).Sxz=bottom(i).Pxz+bottom(i).Dxz+bottom(i).Txz+bottom(i).VDxz+bottom(i).Pixz;
      bottom(i).Syz=bottom(i).Pyz+bottom(i).Dyz+bottom(i).Tyz+bottom(i).VDyz+bottom(i).Piyz;
  end
  
  %With mean convection
  for i=1:length(top)
      top(i).Scxx=top(i).Pxx+top(i).Dxx+top(i).Txx+top(i).VDxx+top(i).Pixx-top(i).Cxx; 
      top(i).Scyy=top(i).Pyy+top(i).Dyy+top(i).Tyy+top(i).VDyy+top(i).Piyy-top(i).Cyy;
      top(i).Sczz=top(i).Pzz+top(i).Dzz+top(i).Tzz+top(i).VDzz+top(i).Pizz-top(i).Czz;
      top(i).Scxy=top(i).Pxy+top(i).Dxy+top(i).Txy+top(i).VDxy+top(i).Pixy-top(i).Cxy; 
      top(i).Scxz=top(i).Pxz+top(i).Dxz+top(i).Txz+top(i).VDxz+top(i).Pixz-top(i).Cxz;
      top(i).Scyz=top(i).Pyz+top(i).Dyz+top(i).Tyz+top(i).VDyz+top(i).Piyz-top(i).Cyz;
  end
  
  for i=1:length(bottom)
      bottom(i).Scxx=bottom(i).Pxx+bottom(i).Dxx+bottom(i).Txx+bottom(i).VDxx+bottom(i).Pixx-bottom(i).Cxx; 
      bottom(i).Scyy=bottom(i).Pyy+bottom(i).Dyy+bottom(i).Tyy+bottom(i).VDyy+bottom(i).Piyy-bottom(i).Cyy;
      bottom(i).Sczz=bottom(i).Pzz+bottom(i).Dzz+bottom(i).Tzz+bottom(i).VDzz+bottom(i).Pizz-bottom(i).Czz;
      bottom(i).Scxy=bottom(i).Pxy+bottom(i).Dxy+bottom(i).Txy+bottom(i).VDxy+bottom(i).Pixy-bottom(i).Cxy; 
      bottom(i).Scxz=bottom(i).Pxz+bottom(i).Dxz+bottom(i).Txz+bottom(i).VDxz+bottom(i).Pixz-bottom(i).Cxz;
      bottom(i).Scyz=bottom(i).Pyz+bottom(i).Dyz+bottom(i).Tyz+bottom(i).VDyz+bottom(i).Piyz-bottom(i).Cyz;
  end
  
  %Skewness tensor. Tensor of Rank 3.
  R3_tensor_tot(1:3,1:3,1:3,1:ln,1:np)=0;
  R3_tensor(1:3,1:3,1:3,1:ln,1:np)=0;
  
  %Form of the tensor.
  % [ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ] 
  % [ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw vww ]
  % [ uuw uvw uww ] [ uvw vvw vww ] [ uww vww www ]
  
  for j=1:np
      for i=1:ln
          R3_tensor_tot(:,:,1,i,j) = [F24(i,j) F28(i,j) F29(i,j); 
                                      F28(i,j) F30(i,j) F34(i,j); 
                                      F29(i,j) F34(i,j) F32(i,j)];
               
          R3_tensor_tot(:,:,2,i,j) = [F28(i,j) F30(i,j) F34(i,j); 
                                      F30(i,j) F25(i,j) F31(i,j); 
                                      F34(i,j) F31(i,j) F33(i,j)];
               
          R3_tensor_tot(:,:,3,i,j) = [F29(i,j) F34(i,j) F32(i,j); 
                                      F34(i,j) F31(i,j) F33(i,j); 
                                      F32(i,j) F33(i,j) F26(i,j)]; 
      end    
  end
  
  for j=1:np
      for i=1:ln
          R3_tensor(:,:,1,i,j)=[(R3_tensor_tot(1,1,1,i,j)-3*U(i,j).*uu(i,j)-U(i,j).*U(i,j).*U(i,j)) ...
                                (R3_tensor_tot(1,2,1,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ... 
                                (R3_tensor_tot(1,3,1,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j));
                                (R3_tensor_tot(1,2,1,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ... 
                                (R3_tensor_tot(2,2,1,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                                (R3_tensor_tot(2,3,1,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j));
                                (R3_tensor_tot(1,3,1,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j)) ...
                                (R3_tensor_tot(2,3,1,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                                (R3_tensor_tot(3,3,1,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j))];
                             
          R3_tensor(:,:,2,i,j)=[(R3_tensor_tot(1,1,2,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ... 
                                (R3_tensor_tot(1,2,2,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                                (R3_tensor_tot(1,3,2,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j));
                                (R3_tensor_tot(1,2,2,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                                (R3_tensor_tot(2,2,2,i,j)-3*V(i,j).*vv(i,j)-V(i,j).*V(i,j).*V(i,j)) ...
                                (R3_tensor_tot(2,3,2,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j));
                                (R3_tensor_tot(1,3,2,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                                (R3_tensor_tot(2,3,2,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j)) ...
                                (R3_tensor_tot(3,3,2,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j))];
          
                             
          R3_tensor(:,:,3,i,j)=[(R3_tensor_tot(1,1,3,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j)) ...
                                (R3_tensor_tot(1,2,3,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                                (R3_tensor_tot(1,3,3,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j));
                                (R3_tensor_tot(1,2,3,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                                (R3_tensor_tot(2,2,3,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j)) ...
                                (R3_tensor_tot(2,3,3,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j));
                                (R3_tensor_tot(1,3,3,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j)) ...
                                (R3_tensor_tot(2,3,3,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j)) ...
                                (R3_tensor_tot(3,3,3,i,j)-3*W(i,j).*ww(i,j)-W(i,j).*W(i,j).*W(i,j))];
      end
  end
  
  for tc=1:topcount
    j=tindicies(tc);alph=top(tc).alpha;
    Ralphap=[cos(alph) sin(alph) 0; ...
        -sin(alph) cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        aabc(1:3,1:3,1:3)=0;
        adef=R3_tensor(:,:,:,i,j);
        
        for aa=1:3
        for bb=1:3    
        for cc=1:3    
        for dd=1:3
        for ee=1:3    
        for ff=1:3
            aabc(aa,bb,cc)=aabc(aa,bb,cc)+Ralphap(aa,dd)*Ralphap(bb,ee) ...
                *Ralphap(cc,ff)*adef(dd,ee,ff);
        end    
        end
        end
        end    
        end
        end
        
        top(tc).uuu(i,1)=aabc(1,1,1);
        top(tc).vvv(i,1)=aabc(2,2,2);
        top(tc).www(i,1)=aabc(3,3,3);
        top(tc).uuv(i,1)=aabc(1,2,1);
        top(tc).uuw(i,1)=aabc(1,3,1);
        top(tc).uvv(i,1)=aabc(2,2,1);
        top(tc).vvw(i,1)=aabc(2,3,2);
        top(tc).uww(i,1)=aabc(3,3,1);
        top(tc).vww(i,1)=aabc(3,3,2);
        top(tc).uvw(i,1)=aabc(2,3,1);
  
    end
  end
  
  for bc=1:botcount
    j=bindicies(bc);alph=bottom(bc).alpha;
    Ralphapp=[cos(alph) sin(alph) 0; ...
        sin(alph) -cos(alph) 0; ...
        0 0 1];
    for i=1:ln
        aabc(1:3,1:3,1:3)=0;
        adef=R3_tensor(:,:,:,i,j);
        
        for aa=1:3
        for bb=1:3    
        for cc=1:3    
        for dd=1:3
        for ee=1:3    
        for ff=1:3
            aabc(aa,bb,cc)=aabc(aa,bb,cc)+Ralphapp(aa,dd)*Ralphapp(bb,ee) ...
                *Ralphapp(cc,ff)*adef(dd,ee,ff);
        end    
        end
        end
        end    
        end
        end
        
        bottom(bc).uuu(i,1)=aabc(1,1,1);
        bottom(bc).vvv(i,1)=aabc(2,2,2);
        bottom(bc).www(i,1)=aabc(3,3,3);
        bottom(bc).uuv(i,1)=aabc(1,2,1);
        bottom(bc).uuw(i,1)=aabc(1,3,1);
        bottom(bc).uvv(i,1)=aabc(2,2,1);
        bottom(bc).vvw(i,1)=aabc(2,3,2);
        bottom(bc).uww(i,1)=aabc(3,3,1);
        bottom(bc).vww(i,1)=aabc(3,3,2);
        bottom(bc).uvw(i,1)=aabc(2,3,1);
    end
    
  end
  
  %for j=1
  %    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
  %        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
  %        0 0 1];
  %        for i=1:ln
  %            aabc(1:3,1:3,1:3)=0;
  %            adef=R3_tensor(:,:,:,i,j);
  %            
  %            for aa=1:3
  %            for bb=1:3    
  %            for cc=1:3    
  %            for dd=1:3
  %            for ee=1:3    
  %            for ff=1:3
  %                aabc(aa,bb,cc)=aabc(aa,bb,cc)+Ralphapp(aa,dd)*Ralphapp(bb,ee) ...
  %                    *Ralphapp(cc,ff)*adef(dd,ee,ff);
  %            end    
  %            end
  %            end
  %            end    
  %            end
  %            end
  %            
  %            bottom(j).uuu(i,1)=aabc(1,1,1);
  %            bottom(j).vvv(i,1)=aabc(2,2,2);
  %            bottom(j).www(i,1)=aabc(3,3,3);
  %            bottom(j).uuv(i,1)=aabc(1,2,1);
  %            bottom(j).uuw(i,1)=aabc(1,3,1);
  %            bottom(j).uvv(i,1)=aabc(2,2,1);
  %            bottom(j).vvw(i,1)=aabc(2,3,2);
  %            bottom(j).uww(i,1)=aabc(3,3,1);
  %            bottom(j).vww(i,1)=aabc(3,3,2);
  %            bottom(j).uvw(i,1)=aabc(2,3,1);
  %
  %        end
  %end
  
  % Some additional terms
  
  if streamwise==1
  
    for j=1:topcount
      shear=nu*rho*gradient(top(j).U,top(j).yn);
      tauw = shear(1);
      ut=sqrt(abs(tauw)/rho);
      yp=ut.*(top(j).yn)./nu;
    
      top(j).tauw=tauw;
      top(j).ut=ut;
      top(j).yp=yp;
      top(j).Up=top(j).U/ut;
      top(j).uup=top(j).uu/ut;
      top(j).vvp=top(j).vv/ut;
      top(j).wwp=top(j).ww/ut;
      top(j).uvp=top(j).uv/ut;
      top(j).vwp=top(j).vw/ut;
  
    end    
      
    for j=1:botcount
      shear=nu*rho*gradient(bottom(j).U,bottom(j).yn);
      tauw = shear(1);
      ut=sqrt(abs(tauw)/rho);
      yp=ut.*(bottom(j).yn)./nu;
    
      bottom(j).tauw=tauw;
      bottom(j).ut=ut;
      bottom(j).yp=yp;
      bottom(j).Up=bottom(j).U/ut;
      bottom(j).uup=bottom(j).uu/ut;
      bottom(j).vvp=bottom(j).vv/ut;
      bottom(j).wwp=bottom(j).ww/ut;
      bottom(j).uvp=bottom(j).uv/ut;
      bottom(j).vwp=bottom(j).vw/ut;
   
    end 
  
  elseif streamwise==3
  
    for j=1:topcount
      shear=nu*rho*gradient(top(j).W,top(j).yn);
      tauw = shear(1);
      ut=sqrt(abs(tauw)/rho);
      yp=ut.*(top(j).yn)./nu;
    
      top(j).tauw=tauw;
      top(j).ut=ut;
      top(j).yp=yp;
      top(j).Up=top(j).U/ut;
      top(j).Wp=top(j).W/ut;
      top(j).uup=top(j).uu/(ut^2);
      top(j).vvp=top(j).vv/(ut^2);
      top(j).wwp=top(j).ww/(ut^2);
      top(j).uvp=top(j).uv/(ut^2);
      top(j).vwp=top(j).vw/(ut^2);
     
    end    
      
    for j=1:botcount
      shear=nu*rho*gradient(bottom(j).W,bottom(j).yn);
      tauw = shear(1);
      ut=sqrt(abs(tauw)/rho);
      yp=ut.*(bottom(j).yn)./nu;
    
      bottom(j).tauw=tauw;
      bottom(j).ut=ut;
      bottom(j).yp=yp;
      bottom(j).Up=bottom(j).U/ut;
      bottom(j).Wp=bottom(j).W/ut;
      bottom(j).uup=bottom(j).uu/ut^2;
      bottom(j).vvp=bottom(j).vv/ut^2;
      bottom(j).wwp=bottom(j).ww/ut^2;
      bottom(j).uvp=bottom(j).uv/ut^2;
      bottom(j).vwp=bottom(j).vw/ut^2;

    end 
  
  end         % streamwise
  
  clearvars -except top bottom svfname fname Rer Domain mnelxyz nel Poly nstat nderiv times timee atime DT nrec Tint npoints wallpts nypts streamwise ssf times timee
  
  save(svfname);
  display([svfname ' saved'])
%  [status result] = system(['mv ' svfname ' re750k_mat20/']);

  clearvars -except ssf

  toc

end  





