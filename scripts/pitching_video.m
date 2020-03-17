%     Creating a pitching video

clear
clc
close all

ifsaab = 1;
ifsave = 1;
outpos = [0.25 0.25 0.6 0.4];

lw = 4.;
cols = lines(1);

if ifsaab

  ifspring = 0;
  axis_x0 = 0.35;
  axis_y0 = 0.034;
  fname = 'ed36F128+14.dat';
  fout  = 'saab_pitch';
  amp = 2.0;

else

  ifspring = 1;
  axis_x0 = 0.25;
  axis_y0 = 0.0;
  fname = 'naca0012.dat';
  fout  = 'naca_pitch';
  amp = 5.0;

end


foil = importdata(fname);

x0 = foil(:,1);
y0 = foil(:,2);

figure(1);
set(gcf,'Units','normalized');
set(gcf,'OuterPosition',outpos);
p0 = plot(x0,y0);
ylim([-0.12 0.12])
xlim([-0.05 1.05])
hold on

if (ifspring)
  nspringpts = 1000;
  theta = linspace(0,7*pi,nspringpts);
  rad   = linspace(0,0.05,nspringpts);
  
  sx = axis_x0 + rad.*cos(theta);
  sy = axis_y0 + rad.*sin(theta);
  
  plot(sx,sy, 'k','LineWidth',3.0)
else
  plot(axis_x0,axis_y0, 'ok','LineWidth',2.5,'MarkerSize',12)

end

%return
set(gca,'Visible','off')

fps = 20;
nsecs = 30;
ntheta = fps*nsecs;
theta2 = linspace(0,10*pi,ntheta);
alpha  = amp*pi/180*sin(theta2);

for i=1:ntheta

   delete(p0)

   R = [cos(alpha(i)) sin(alpha(i)); ...
        -sin(alpha(i)) cos(alpha(i))];
   coords = R*[(x0-axis_x0)'; (y0-axis_y0)'];
   x2 = coords(1,:)' + axis_x0;
   y2 = coords(2,:)' + axis_y0;
   
   p0 = plot(x2,y2, 'LineWidth', lw, 'Color', cols(1,:));

   fname2 = sprintf('%s%5.5d.png',fout,i);

   pause(0.01)
   
   if (ifsave)
     SaveFig(gcf,fname2,'plots/',1)
   end  
end








