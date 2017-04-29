% Splinning foil data

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

dx=diff(x);
dy=diff(y);
ds=sqrt(dx.^2 + dy.^2);
s=cumsum(ds);
s = [0; s];

st=s(1:ind);
sb=s(ind:end);

sphase=pi*linspace(0,1,1000);
sfine1=max(st)/2*(1 + cos(sphase));
sfine2=max(sfine1)+(max(sb)-min(sb))/2*(1+cos(sphase));

sfineall = [sfine1(end:-1:2) sfine2(end:-1:1)];


xintp_c = interp1(s,x,sfineall,'cubic');
yintp_c = interp1(s,y,sfineall,'cubic');

xintp_c = (xintp_c-min(xintp_c))/(max(xintp_c) - min(xintp_c));
yintp_c = yintp_c/(max(xintp_c) - min(xintp_c));

xintp_s = interp1(s,x,sfineall,'spline');
yintp_s = interp1(s,y,sfineall,'spline');

xintp_s = (xintp_s-min(xintp_s))/(max(xintp_s) - min(xintp_s));
yintp_s = yintp_s/(max(xintp_s) - min(xintp_s));

xintp_cs = csaps(s,x,0.99999999,sfineall);
yintp_cs = csaps(s,y,0.99999999,sfineall);

xintp_cs = (xintp_cs-min(xintp_cs))/(max(xintp_cs) - min(xintp_cs));
yintp_cs = yintp_cs/(max(xintp_cs) - min(xintp_cs));

% barryx = barylag([s, x],sfineall');
% barryy = barylag([s, y],sfineall');


figure(1)
plot(xintp_c,yintp_c);hold on
plot(xintp_s,yintp_s, 'r')
plot(xintp_cs,yintp_cs, 'k')
%plot(barryx,barryy, 'k')

plot(x,y, '-om', 'LineWidth', 1.5)

kc=CalcCurvature(xintp_c,yintp_c);
ks=CalcCurvature(xintp_s,yintp_s);
kcs=CalcCurvature(xintp_cs,yintp_cs);

figure(2)
plot(xintp_c,kc, '-'); hold on
plot(xintp_s,ks, '-r')
plot(xintp_cs,kcs, '-k')
ylim([0 100])


%% Find residual @ collocated points.
for i=1:length(x)
   x0=x(i);
   y0=y(i);
   dmin=100;
   for k=1:length(xintp_cs);
     dl=sqrt((xintp_cs(k)-x0)^2 + (yintp_cs(k)-y0)^2);
     if dl<dmin
       dmin=dl;
     end
   end    
   resid(i)=dmin;
end

figure(3)
plot(x,resid, '.')   
      
%% Write out data
l1=length(xintp_cs);
fid = fopen('saab_foil2.dat','w');
fprintf(fid,'%1d\t 2\n', l1);
for i=1:l1
  fprintf(fid,'%10.7f\t %10.7f\t 0.00\n',xintp_cs(i),yintp_cs(i));
end
fclose(fid)



