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

%Initial estimation for number of elements
%Start with Chebyshev grid
%Nel=20;         
%19

%Length of computational domain Lx (homogeneous direction)
% Lz=25;
Lz=2*pi;
Lx=pi;

%Criterion for maximum Delta x+
%Standard criterion is Delta x+_max <10; for LES Delta x+_max < 15
dzMaxElPlus=20;                             
dxMaxElPlus=9;                             


firstpt = 0.6;           % Y+ of the first point
%Another criterion is to keep a certain number of points below a particular
%reference close to the wall
%Standard criterion is 4 points (including the wall) below y+=1
%Reference location close to the wall
drWallElPlus=6.0;        % Right now arbitrary for LES 

%Then there is another criterion: 10 points (including the wall) below
%y+=10 (this was original Spalart's criterion. In his retraction, he
%suggests that the 10th point should be at y+=6 or even below).

%Number of points to be kept below reference location
n_wall=6;                % Right now arbitrary for LES

%Last criterion is the maximum allowed Delta y+
%Standard value is Delta y+_max<5
%This value should be modified to obtain the right behavior (with 5, we
%will get more than 5, so it is necessary to use 4.5 or similar)
dthMaxElPlus=9;             % <9 forLES 

%We consider the first part of the mesh to follow a Chebyshev distribution,
%then a uniform spacing is considered towards the center of the duct. 
%Here we chosse the number of elements following Chebyshev distribution, a
%standard value being 5, but it is necessary to check final distribution
x_cheb_reliable=5;

%GRID IN STREAMWISE DIRECTION
%x represens GLL Nodes
[x]=lglnodes(N);  
dz_gll=abs(diff(x));
dz_max=max(dz_gll);
dz_min = min(dz_gll);

dz_max_el=dzMaxElPlus/Re_tau;
len_el_z=dz_max_el/dz_max*2;

dx_max_el=dxMaxElPlus/Re_tau;
len_el_x=dx_max_el/dz_max*2;

%Number of stages needed to maintain the streamwise grid criterion 
nelz=Lz/len_el_z;
nelx=Lx/len_el_x;

%GRID IN WALL-NORMAL DIRECTION
%dr=abs(x/2-0.5);
dy_gll = x;
dr = 1 + dy_gll(end:-1:1);

%Keep certain points below reference
y_ref_n=dr(n_wall);

y_el_n=drWallElPlus/Re_tau;

el_min=2*y_el_n/y_ref_n;

%Check y+ for first point.
y1plus = el_min/2*dz_min*Re_tau;

if y1plus>firstpt           % we want first point at y+<0.95
  dr_min    = min(diff(dr));
  dr_min_el = firstpt/Re_tau;
  el_min    = 2*dr_min_el/dr_min;         % max element size at the wall
end

dy_ref=abs(diff(x));
dy_max_ref=max(dy_ref);
dy_max_el=dthMaxElPlus/Re_tau;
el_max=2*dy_max_el/dy_max_ref;           % max element size based on deltay+max

%We need to modify the Nel value to obtain a mesh satisfying all the
%constrains simultaneously

Nel = 2;
ifbreak=0;

while (1) 

[xel, DM] = chebdif(Nel, 1);

%Distance between Chebychev nodes in element [-1,1]
dx_cheb_ref=abs(diff(xel)); 
dx_cheb_ratio = dx_cheb_ref(2:end)./dx_cheb_ref(1:end-1);
dx_cheb_min=min(dx_cheb_ref);
[dx_cheb_max maxind]=max(dx_cheb_ref);
rtio = dx_cheb_max/dx_cheb_min;

dx_cheb_max_el = rtio*el_min;

%scaled_d_max = len_el_r*d_max_line/d_min_line
el2 = 2*el_min/dx_cheb_min;
dx_cheb_el = dx_cheb_ref*el2/2;

tlength = el2/2*sum(dx_cheb_ref(1:maxind));

if ifbreak
     break;
end


if dx_cheb_max_el>el_max
     Nel=Nel-1;
     ifbreak = 1; 
     continue;
end

if tlength > 1
     Nel=Nel-1;
     ifbreak = 1;
     continue;
end

Nel=Nel+1;

end

el2 = 2*el_min/dx_cheb_min;
dx_cheb_el = dx_cheb_ref*el2/2;

%% Require minimum number of Chebeshev elements: 15 for now
nelmin = 7;

if Nel<nelmin
  Nel = nelmin;
 
  [xel, DM] = chebdif(Nel, 1);

  %Distance between Chebychev nodes in element [-1,1]
  dx_cheb_ref=abs(diff(xel)); 
  dx_cheb_ratio = dx_cheb_ref(2:end)./dx_cheb_ref(1:end-1);
  dx_cheb_min=min(dx_cheb_ref);
  [dx_cheb_max maxind]=max(dx_cheb_ref);
  rtio = dx_cheb_max/dx_cheb_min;

  dx_cheb_min_el = el_max/rtio;

  el2 = 2*el_max/dx_cheb_max;
  dx_cheb_el = dx_cheb_ref*el2/2;

  tlength = el2/2*sum(dx_cheb_ref(1:maxind));
end 


%Combine Chebyshev and uniform distribution
nodes(1) = 1;
for i = 2:maxind+1
  nodes(i) = nodes(i-1)-dx_cheb_el(i-1);
end

%break

length_equi=1-tlength;
twice_len = length_equi*2;
nel_equi = twice_len/el_max;
needed_equi=ceil(nel_equi);
% new_ratio=length_tobe_covered_equi/needed_equi;
new_delta = twice_len/needed_equi;

allnodes = nodes;
l1=length(allnodes);
for i = 1:needed_equi
  allnodes(l1+i) = allnodes(l1+i-1)-new_delta;
end

%break
allnodes = [allnodes -nodes(end-1:-1:1)];
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

figure(2);
plot(mesh_final,'ob'); hold on;
plot(moser(2).yp,'or'); hold on;
% plot(zz,'og'); hold on;
legend('Current','Moser','Location','NorthWest');
N_total=number_of_elements^2*nel_z*(N+1)^3;

figure(3);
plot(diff(mesh_final), 'k')

%% save variables for box file
nely = number_of_elements;
Ly = 2;
nelx = ceil(nelx);
nelz = ceil(nelz);
ypts = allnodes(end:-1:1);

%close all

clearvars -except nelx nelz nely Lx Lz Ly ypts
save boxpar
