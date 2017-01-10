clc 
close all
warning off all
clear all
format long

load 'channel_moser.mat'
% load 'test_18_330.mat'
% data = importdata('chan395.means');

%INPUT PARAMETERS
%Reynolds number based on the friction velocity
Re_tau=395;     

%Polynomial order. (Number of points on the element -1)
N=11;            

%Length of computational domain Lx (homogeneous direction)
%Lz=25;
Lz=2*pi;
Lx=pi;

%Criterion for maximum Delta x+
%Standard criterion is Delta x+_max <10; for LES Delta x+_max < 15
avgdz = 1;
dzMaxElPlus=20;                             

avgdx = 1;
dxMaxElPlus=7.1;                             

firstpt = 0.64;           % Y+ of the first point
%Another criterion is to keep a certain number of points below a particular
%reference close to the wall
%Standard criterion is 4 points (including the wall) below y+=1
%Reference location close to the wall
drWallElPlus=6.0;        % Right now arbitrary for LES 

%Then there is another criterion: 10 points (including the wall) below
%y+=10 (this was original Spalart's criterion. In his retraction, he
%suggests that the 10th point should be at y+=6 or even below).
%Number of points to be kept below reference location
n_wall=1;                % Right now arbitrary for LES

%Last criterion is the maximum allowed Delta y+
%Standard value is Delta y+_max<5
%This value should be modified to obtain the right behavior (with 5, we
%will get more than 5, so it is necessary to use 4.5 or similar)
dthMaxElPlus=12;             % < arbitrary for les right now 

%GRID IN STREAMWISE DIRECTION
%x represens GLL Nodes
[x]=lglnodes(N);  
dz_gll=abs(diff(x));
dz_max=max(dz_gll);
dz_min = min(dz_gll);
dz_mean = mean(dz_gll);

dz_max_el = dzMaxElPlus/Re_tau;
len_el_z=dz_max_el/dz_max*2;
nelz=Lz/len_el_z;

if avgdz
     dz_mean_el = dzMaxElPlus/Re_tau;
     len_el_z=dz_max_el/dz_mean*2;
     nelz = Lz/len_el_z;
     nelz=ceil(nelz);
     len_el_z = Lz/nelz;
end
maxdzplus = len_el_z/2*dz_max*Re_tau;
meandzplus = len_el_z/2*mean(dz_gll)*Re_tau;

%Number of stages needed to maintain the streamwise grid criterion 
[x]=lglnodes(N);  
dx_gll=abs(diff(x));
dx_max=max(dx_gll);
dx_min = min(dx_gll);
dx_mean = mean(dx_gll);

dx_max_el=dxMaxElPlus/Re_tau;
len_el_x=dx_max_el/dx_max*2;
if avgdx
     dx_mean_el = dxMaxElPlus/Re_tau;
     len_el_x=dx_max_el/dx_mean*2;
     nelx = Lx/len_el_x;
     nelx=ceil(nelx);
     len_el_x = Lx/nelx;
end
nelx=Lx/len_el_x;
maxdxplus = len_el_x/2*dx_max*Re_tau;
meandxplus = len_el_x/2*mean(dx_gll)*Re_tau;

%GRID IN WALL-NORMAL DIRECTION


dy_gll = x;
dr = 1 + dy_gll(end:-1:1);

dy_ref=dz_gll;
dy_ref_min = min(dy_ref);
dy_ref_max = max(dy_ref);
dy_max_el=dthMaxElPlus/Re_tau;

y1st =firstpt/Re_tau;

el1st = 2*y1st/dy_ref_min;
el_max = 2*dy_max_el/dy_ref_max;           % max element size based on deltay+max

% Initial ratios r1,r2
r1 = 1.04;
r2 = 1.06;
max_r = 1.6;
ini_lin = 0;
ini_lin_plus = 50;

cond = 1;
nely = 1;
ypts = [0 el1st];
el_ratio = [1];
tlength = max(ypts);
count = 0;

prev_size = el1st;

while (cond)
count = count+1;

ratio = r1*(r2^count);
if ratio > max_r 
     ratio = max_r;                          % Maximum ratio  
end

if (ini_lin)                                   % Restrictions for the first few elements. (iftrue)
     tlen_plus = tlength*Re_tau;
     if tlen_plus<ini_lin_plus;
          ratio=r1;
          count=count-1;
     end
end

nxt_el_size = prev_size*ratio;

if nxt_el_size>el_max
     break;
end
tl = tlength + nxt_el_size;
if tl>1
     break;
end

tlength = tl;
ypts(end+1) = ypts(end) + nxt_el_size;
el_ratio(end+1) = nxt_el_size/prev_size;
prev_size = nxt_el_size;

end

tleft = 2*(1-tlength);
nlin = tleft/el_max;

if nlin>0.8 && nlin<=1.0
     fin_size = tleft;
     linpts = tlength + fin_size;

     ypts = ypts -1;
     linpts = linpts-1;
     nodes = [ypts linpts -ypts(end-1:-1:1)];
    

elseif nlin<2
     ypts(end) = [];
     tlength = ypts(end);
     tleft = 2*(1-tlength);
     nlin = tleft/el_max;
     nlin = ceil(nlin);

     for i = 1:nlin
          linpts(i) = ypts(end) + i*tleft/nlin;
     end

     ypts = ypts-1;
     linpts = linpts -1;
     nodes = [ypts linpts -ypts(end-1:-1:1)];
     
else

     nlin = ceil(nlin);

     for i = 1:nlin
          linpts(i) = ypts(end) + i*tleft/nlin;
     end

     ypts = ypts-1;
     linpts = linpts -1;
     nodes = [ypts linpts -ypts(end-1:-1:1)];

end

Nel = size(nodes);

%break
allnodes = nodes 
nodes=allnodes';
%%

%x_mid=zeros(needed_equi+1,1);
%x_mid(1)=length_tobe_covered_equi;
%for ii=1:needed_equi
%    x_mid(ii+1)=x_mid(ii)-new_ratio;
%end
%nodes=[xel(1:x_cheb_reliable);x_mid;-flipud(x_mid);xel(end-x_cheb_reliable+1:end)];

nodes_r=unique(nodes)';
number_of_nodes=length(nodes_r);
number_of_elements=number_of_nodes-1;

%nodes_final=zeros(number_of_nodes,1);

%count=0;
%for i=1:length(nodes_r)
%    if abs(nodes_r(i))>1e-6
%        count=count+1;
%        nodes_final(count)=nodes_r(i);
%    else
%        if nodes_r(i)>0
%            count=count+1;
%            nodes_final(count)=0;
%        end
%    end
%end

nodes_final=nodes_r;
number_of_nodes;
number_of_elements;

gll_nodes=lglnodes(N)/2+0.5;

delta_yw_plus=abs((nodes_final(1)-nodes_final(2)))*Re_tau*gll_nodes(end-1);
disp(['delta y+ at the wall =', num2str(delta_yw_plus)]);
% largest y+
max_el_size = max(abs(diff(nodes_final)));
delta_yc_max = max_el_size*Re_tau*max(abs(diff(gll_nodes)));
disp(['Max delta y+ =', num2str(delta_yc_max)]);
%delta_yc_plus=(nodes(round(length(nodes)/2-1)))*Re_tau*(gll_nodes(end/2)-gll_nodes(end/2+1))


yplus(1)=0;
yplus(2)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-1); 
yplus(3)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-2);
yplus(4)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-3);
yplus(5)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-4);
yplus(6)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-5);
yplus(7)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-6);
yplus(8)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-7);
yplus(9)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-8);
yplus(10)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-9);
yplus(11)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-10);
yplus(12)=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-11);
disp('Node no                 Delta y+')
disp([[1:length(yplus)]' yplus'])

nel_z=nelz;
delta_z_avg=Lz/((N+1)*nel_z)*Re_tau;
delta_z_max=delta_z_avg*(gll_nodes(end/2)-gll_nodes(end/2+1))/(1/(N+1));
delta_z_min=delta_z_avg*gll_nodes(end-1)/(1/(N+1));

nel_x=nelx;
delta_x_avg=Lx/((N+1)*nel_x)*Re_tau;
delta_x_max=delta_x_avg*(gll_nodes(end/2)-gll_nodes(end/2+1))/(1/(N+1));
delta_x_min=delta_x_avg*gll_nodes(end-1)/(1/(N+1));

for i=1:number_of_elements
    dy_final(i)=nodes_final(i+1)-nodes_final(i);
end

figure(1);
plot(nodes_final,'o-b'); hold on;
plot([0 dy_final],'o-r'); hold on;

y_last=0;
for i=1:number_of_elements/2
    for j=1:N+1
        mesh_final((i-1)*(N+1)+j,1)=gll_nodes(end-j+1)*dy_final(i)*Re_tau+y_last;
        if j==N+1 
            y_last=gll_nodes(end-j+1)*dy_final(i)*Re_tau+y_last;
        end
    end
end

%mesh_final = unique(mesh_final); % Consecutive elements share points
figure(2);
plot(mesh_final,'ob'); hold on;
plot(moser(2).yp,'or'); hold on;
% plot(zz,'og'); hold on;
legend('Current','Moser','Location','NorthWest');
N_total=number_of_elements^2*nel_z*(N+1)^3;

figure(3);
plot(diff(mesh_final), 'k')

el_sizes = abs(diff(nodes_final));
el_ratio = [1];
for j = 1:length(el_sizes)-1
     el_ratio(j+1) = el_sizes(j+1)/el_sizes(j);
end

figure(4)
plot(el_ratio, '-*', 'MarkerSize', 16)
title('Ratios of consecutive sizes')

%% save variables for box file
nely = number_of_elements;
Ly = 2;
nelx = ceil(nelx);
nelz = ceil(nelz);
ypts = allnodes;

maxdzplus
meandzplus
maxdxplus
meandxplus

%close all

clearvars -except nelx nelz nely Lx Lz Ly ypts
nelx
nely
nelz
save boxpar
