clear
clc
close all

lafs=20;

%---------------------------------------- 
snormals = importdata('surf_normals.67');

x_imp11 = snormals.data(:,1);
[x_imp11 I] = sort(x_imp11);
y_imp11 = snormals.data(I,2);

snx11 = -snormals.data(I,4);
sny11 = -snormals.data(I,5);

ind = sny11>0;
snx_top11 = snx11(find(ind));
sny_top11 = sny11(find(ind));
xt_imp11  = x_imp11(find(ind));
yt_imp11  = y_imp11(find(ind));

snx_bot11 = snx11(find(~ind));
sny_bot11 = sny11(find(~ind));
xb_imp11  = x_imp11(find(~ind));
yb_imp11  = y_imp11(find(~ind));

for i=1:length(snx_top11)
  stx_top11(i) = sny_top11(i);
  sty_top11(i) = -snx_top11(i);
end

for i=1:length(snx_bot11)
  stx_bot11(i) = -sny_bot11(i);
  sty_bot11(i) = snx_bot11(i);
end

%% Rotate coordinates

% Rotate imported values according to simulation time   
dtheta = -6.7*pi/180;
axis_x0=0.35;
axis_y0=0.034;

% Top side
xnew = xt_imp11;
ynew = yt_imp11;

% positive clockwise              
rot = [cos(dtheta) sin(dtheta); ...
       -sin(dtheta) cos(dtheta)];

coords = rot*[transpose(xnew)-axis_x0; transpose(ynew)-axis_y0];
xrtnew = coords(1,:) + axis_x0;
yrtnew = coords(2,:) + axis_y0;

% Bottom side
xnew = xb_imp11;
ynew = yb_imp11;

% positive clockwise              
rot = [cos(dtheta) sin(dtheta); ...
       -sin(dtheta) cos(dtheta)];

coords = rot*[transpose(xnew)-axis_x0; transpose(ynew)-axis_y0];
xrbnew = coords(1,:) + axis_x0;
yrbnew = coords(2,:) + axis_y0;


%% end of rotation   

plot(xrtnew,yrtnew, 'LineWidth', 3); hold on
plot(xrbnew,yrbnew, 'LineWidth', 3);
xlim([0 1])
ylim([-0.1 0.15])
pbaspect([4 1 1])

xlabel('$x$', 'FontSize', lafs)
ylabel('$y$', 'FontSize', lafs)
SaveFig(gcf, 'foil.eps', 'plots/', 1)







