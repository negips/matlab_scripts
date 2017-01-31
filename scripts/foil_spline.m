% splining foil data

clear
clc
close all

foil = importdata('ed36F128+14.dat');
x = foil.data(:,1);
y = foil.data(:,2);

[xmin ind] = min(x);
topx = x(1:ind);
topy = y(1:ind);

botx = x(ind:end);
boty = y(ind:end);

%topx = [topx; botx(2)];
%topy = [topy; boty(2)];
le_ns = 3;                    % Special treatment for leading edge points

plot(topx,topy, '-sb')
hold on
plot(botx,boty, '-sr')

nseg = length(topx)-1;

topfinex = [];
topfiney_cu = [];
topfiney_sp = [];

for i=1:nseg-le_ns                % Skip the last point
  x0=topx(i);
  x1= topx(i+1);

  if (x1<0.2)
    nreso=20;
  else
    nreso=10;
  end                
  
  xintp = linspace(x0,x1,nreso);
  xintp = xintp(1:end-1);                              % Remove last point
  yintp_cu = interp1(topx,topy,xintp,'cubic');
  yintp_sp = spline(topx,topy,xintp);

  topfinex = [topfinex xintp];
  topfiney_cu = [topfiney_cu yintp_cu];    
  topfiney_sp = [topfiney_sp yintp_sp];    

end

plot(topfinex,topfiney_cu, '*k')
%plot(topfinex,topfiney_sp, 'dm')


nseg = length(botx)-1;

botfinex = [];
botfiney_cu = [];
botfiney_sp = [];

for i=1+le_ns:nseg                  % Skip the first point
  x0=botx(i);
  x1= botx(i+1);

  if (x1<0.2)
    nreso=20;
  else
    nreso=10;
  end                
 
  xintp = linspace(x0,x1,nreso);

  if i~=nseg
    xintp = xintp(1:end-1);         % Remove last point
  end                    
  yintp_cu = interp1(botx,boty,xintp,'cubic');
  yintp_sp = spline(botx,boty,xintp);    

  botfinex = [botfinex xintp];
  botfiney_cu = [botfiney_cu yintp_cu];    
  botfiney_sp = [botfiney_sp yintp_sp];    

end

plot(botfinex,botfiney_cu, '*k')
%plot(botfinex,botfiney_sp, 'dm')


%% Leading edge

xvec = [topx(end-le_ns:end) ; botx(2:le_ns+1)];
yvec = [topy(end-le_ns:end) ; boty(2:le_ns+1)];
% For these points, do interpolation on x instead of y

nseg = length(yvec)-1;

lefiney = [];
lefinex_cu = [];
lefinex_sp = [];

nreso = 20;
for i=1:nseg
  y0 = yvec(i);
  y1 = yvec(i+1);    

  yintp = linspace(y0,y1,nreso);
  yintp = yintp(1:end-1);                 % Remove last point
  xintp_cu = interp1(yvec,xvec,yintp,'cubic');
  xintp_sp = spline(yvec,xvec,xintp);

  lefiney = [lefiney yintp];
  lefinex_cu = [lefinex_cu xintp_cu];    
  lefinex_sp = [lefinex_sp xintp_sp];

end 

plot(lefinex_cu,lefiney, 'sg')


xfineall = [topfinex lefinex_cu botfinex];
yfineall = [topfiney_cu lefiney botfiney_cu];

%[xq,IA,IC] = unique(xfineall);

l1=length(xfineall);
fid = fopen('saab_foil.dat','w');
fprintf(fid,'%1d\t 2\n', l1);
for i=1:l1
  fprintf(fid,'%10.7f\t %10.7f\t 0.00\n',xfineall(i),yfineall(i));
end
fclose(fid)


