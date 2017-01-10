%----------------------------------------------------------------------
%write initial conditions from RANS solution
%----------------------------------------------------------------------
disp('Write initial conditions')
Nelm2D = length(EL)/nelz;
Nelm = length(EL);
%gll2D = zeros(length(EL)*(N+1)^2,2);
gll2D = zeros(Nelm2D*(N+1)^2,2);
for i = 1:Nelm2D
    gll2D((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL(i).GLL(1:(N+1)^2,1:2);
end

disp('Interpolatin ...')
%interpolate RANS solution onto gll points
uinit = griddata(xr,yr,ur,gll2D(:,1),gll2D(:,2));
vinit = griddata(xr,yr,vr,gll2D(:,1),gll2D(:,2));
winit = griddata(xr,yr,wr,gll2D(:,1),gll2D(:,2));
pinit = griddata(xr,yr,pr,gll2D(:,1),gll2D(:,2));

%%% Smooth fields to remove "interpolation edginess"
%swidth = 29;
%uinit = smooth(uinit,swidth);
%vinit = smooth(vinit,swidth);
%winit = smooth(winit,swidth);
%pinit = smooth(pinit,swidth);

disp('Copying ...')
data = zeros(length(EL),(N+1)^3,7);
for i = 1:length(EL)
    if (mod(i,10000)==0)
        disp([ num2str(length(EL)-i) ' elements left'  ]);
    end
       
    for j = 1:Ngp
        if (j==1 && (i <= Nelm2D )) || ( j ==Ngp &&  i>= (Nelm - Nelm2D) )
          index2D = (i-1)*(N+1)^2+1-floor((i-1)/Nelm2D)*Nelm2D*(N+1)^2:i*(N+1)^2-floor((i-1)/Nelm2D)*Nelm2D*(N+1)^2;
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,1) = uinit(index2D);
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,2) = vinit(index2D);
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,3) = winit(index2D);
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,4) = pinit(index2D);
        else 
            
        index2D = (i-1)*(N+1)^2+1-floor((i-1)/Nelm2D)*Nelm2D*(N+1)^2:i*(N+1)^2-floor((i-1)/Nelm2D)*Nelm2D*(N+1)^2;
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,1) = uinit(index2D);
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,2) = vinit(index2D);
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,3) = winit(index2D);
        EL(i).Vinit((j-1)*(N+1)^2+1:j*(N+1)^2,4) = pinit(index2D);
        end
    end
    
    %write data-file to be written out
    data(i,:,1:3) = EL(i).GLL;
    data(i,:,4:7) = EL(i).Vinit;
end

N1 = N+1;

%write nekton restart file
%Header: 2D field
%hdr = '#std 8        1                          0.00000000000E+00         0      0      1 XUP                                              ';
hdr = '#std 4        1                          0.00000000000E+00         0      0      1 XUP                                              ';
%hdr = '#std 4        1                        0.0000000000000E+00         0      0      1 XUP                                              ';

% N:
Nst = num2str(N1);
hdr(9-length(Nst)+1:9)   = Nst;
hdr(12-length(Nst)+1:12) = Nst;
hdr(15-length(Nst)+1:15) = Nst;
% nel:
Nelmst = num2str(Nelm);
hdr(26-length(Nelmst)+1:26) = Nelmst;
hdr(37-length(Nelmst)+1:37) = Nelmst;

%%% tag
tag = 6.54321;
disp('Writing ...')
%stsw = schreib(data,Nelm,hdr,tag,1,3,'le',[path casename '.IC'],'xup');
tag = 6.54321;

stsw = schreib_sp(data,Nelm,hdr,tag,1,3,'le',[path casename '.IC'],'xup');

%clear data;

%clear xb yb ub vb wb;
%clear xr yr ur vr wr;
