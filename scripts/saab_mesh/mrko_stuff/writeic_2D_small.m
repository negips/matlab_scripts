%----------------------------------------------------------------------
%write initial conditions from nek
%----------------------------------------------------------------------
disp('Write initial conditions')

%% import and reshape
gen_import_gllsmall

gll2Dsmall = zeros(length(ELsmall)*(N+1)^2,2);
for i = 1:length(ELsmall)
    gll2Dsmall((i-1)*(N+1)^2+1:i*(N+1)^2,:) = ELsmall(i).GLL(1:(N+1)^2,1:2);
end
%% import baseflow data
gen_import_uvp

%%
% tri = delaunay(xr,yr);
% h = trisurf(tri, xr, yr, pr);
% view(2)
% shading interp
% colorbar EastOutside
% axis equal

%% interpolate baseflow solution on smaller grid
method = 'natural';
uinit = griddata(gll2Dbase(:,1),gll2Dbase(:,2),vel2Dbase(:,1),gll2Dsmall(:,1),gll2Dsmall(:,2),method);
vinit = griddata(gll2Dbase(:,1),gll2Dbase(:,2),vel2Dbase(:,2),gll2Dsmall(:,1),gll2Dsmall(:,2),method);
winit = griddata(gll2Dbase(:,1),gll2Dbase(:,2),vel2Dbase(:,3),gll2Dsmall(:,1),gll2Dsmall(:,2),method);
pinit = griddata(gll2Dbase(:,1),gll2Dbase(:,2),p2Dbase,gll2Dsmall(:,1),gll2Dsmall(:,2),method);
sum(isnan(uinit))
sum(isnan(vinit))
sum(isnan(winit))
sum(isnan(pinit))
%% fit interpolated data in element structure for writing nek input file
data = zeros(length(ELsmall),(N+1)^2,6);
for i = 1:length(ELsmall)
    for j = 1:1 %2D fields
        ELsmall(i).VEL((j-1)*(N+1)^2+1:j*(N+1)^2,1) = uinit((i-1)*(N+1)^2+1:i*(N+1)^2);
        ELsmall(i).VEL((j-1)*(N+1)^2+1:j*(N+1)^2,2) = vinit((i-1)*(N+1)^2+1:i*(N+1)^2);
        ELsmall(i).VEL((j-1)*(N+1)^2+1:j*(N+1)^2,3) = winit((i-1)*(N+1)^2+1:i*(N+1)^2);
        ELsmall(i).P((j-1)*(N+1)^2+1:j*(N+1)^2,1) = pinit((i-1)*(N+1)^2+1:i*(N+1)^2);
    end
%     if(i < box(1).nelx*(box(1).nely-1))
%         EL(i).Vinit(:,1) = .50;
%         EL(i).Vinit(:,2) = 0.0;
%     end
    if(i < box(1).nelx+1)
        ELsmall(i).VEL(1:N+1,1) = 0.0;
        ELsmall(i).VEL(1:N+1,2) = 0.0;
        ELsmall(i).VEL(1:N+1,3) = 0.0;
    end
    %write data-file to be written out
%     data(i,:,1:2) = ELsmall(i).GLL;
%     data(i,:,3:5) = ELsmall(i).VEL;
%     data(i,:,6) = ELsmall(i).P;
    
    data(i,:,1) = ELsmall(i).GLL(:,1);
    data(i,:,2) = ELsmall(i).GLL(:,2);
    data(i,:,3) = ELsmall(i).VEL(:,1);
    data(i,:,4) = ELsmall(i).VEL(:,2);
    data(i,:,5) = ELsmall(i).P(:,1);
    data(i,:,6) = ELsmall(i).VEL(:,3);
end

% %% write nek input file
% N1 = N+1;
% %Header: 2D field
% hdr = '#std 8        1                         0.000000000000E+00        0      0      1 XUPT                                              ';
% % N:
% Nst = num2str(N1);
% Nstz = num2str(1);
% hdr(9-length(Nst)+1:9)   = Nst;
% hdr(12-length(Nst)+1:12) = Nst;
% hdr(15-length(Nstz)+1:15) = Nstz;
% % nel:
% Nelmst = num2str(Nelm);
% hdr(26-length(Nelmst)+1:26) = Nelmst;
% hdr(37-length(Nelmst)+1:37) = Nelmst;
% 
% %%% tag
% tag = 6.54321;
% stsw = schreib(data,Nelm,hdr,tag,1,2,'le',[pathsmall casename '.IC'],'xupt');
% 
% clear data;
% disp('Write boundary conditions')

%% check with plots
[xx, yy, uu, vv, pp] = rearrange_2d(ELsmall,box(1).nelx,box(1).nely,N1);

%%
figure
plot(pp(:,1),yy(:,1),'b*-',pp(:,end),yy(:,end),'r*-');grid on
legend('lower out flow','upper out flow','location','northwest')
xlabel('p'),ylabel('y')
[maxu_lo,indmaxu_lo] = max(pp(:,1));
[maxu_up,indmaxu_up] = max(pp(:,end));

%%
figure
pcolor(xx,yy,uu);shading interp;colorbar;axis equal
xlabel('x'),ylabel('y'),title('u')

%%
figure
pcolor(xx,yy,pp);shading interp;colorbar;axis equal
xlabel('x'),ylabel('y'),title('p')








