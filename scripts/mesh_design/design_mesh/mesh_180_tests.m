clc 
close all
warning off all
clear all
format long

load 'channel_moser.mat'

%INPUT PARAMETERS
%Reynolds number based on the friction velocity
Re_tau=180;     

%Polynomial order. (Number of points on the element -1)
N=11;            

%Initial estimation for number of elements
%Start with Chebyshev grid
Nel=14;         
%19

%Length of computational domain Lx (homogeneous direction)
Lz=4*pi;

%Criterion for maximum Delta x+
%Standard criterion is Delta x+_max <10
dzMaxElPlus=10; 

%Another criterion is to keep a certain number of points below a particular
%reference close to the wall
%Standard criterion is 4 points (including the wall) below y+=1
%Reference location close to the wall
drWallElPlus=1.0; 

%Then there is another criterion: 10 points (including the wall) below
%y+=10 (this was original Spalart's criterion. In his retraction, he
%suggests that the 10th point should be at y+=6 or even below).

%Number of points to be kept below reference location
n_wall=4;

%Last criterion is the maximum allowed Delta y+
%Standard value is Delta y+_max<5
%This value should be modified to obtain the right behavior (with 5, we
%will get more than 5, so it is necessary to use 4.5 or similar)
dthMaxElPlus=4.2; 

%We consider the first part of the mesh to follow a Chebyshev distribution,
%then a uniform spacing is considered towards the center of the duct. 
%Here we chosse the number of elements following Chebyshev distribution, a
%standard value being 5, but it is necessary to check final distribution
x_cheb_reliable=4;

%GRID IN STREAMWISE DIRECTION
%x represens GLL Nodes
[x]=lglnodes(N);  
dz=abs(x(2:end)-x(1:end-1));
dz_max=max(dz)/2;

dz_max_el=dzMaxElPlus/Re_tau;
len_el_z=dz_max_el/dz_max;

%Number of stages needed to maintain the streamwise grid criterion 
nelz=Lz/len_el_z;

%GRID IN WALL-NORMAL DIRECTION
dr=abs(x/2-0.5);

%Keep 4 points below y+=1
dr_min=dr(n_wall);  
dr_min_el=drWallElPlus/Re_tau;
len_el_r=dr_min_el/dr_min;

dth=abs(x(2:end)-x(1:end-1));
dth_max=max(dth)/2;
dth_max_el=dthMaxElPlus/Re_tau;
len_el_th=dth_max_el/dth_max;

%We need to modify the Nel value to obtain a mesh satisfying all the
%constrains simultaneously
[xel, DM] = chebdif(Nel, 1);

%Distance between Chebychev nodes in element [-1,1]
d_line=abs(xel(2:end)-xel(1:end-1)); 
d_min_line=min(d_line);
d_max_line=max(d_line);

%Combine Chebyshev and uniform distribution
length_tobe_covered_equi=xel(x_cheb_reliable);
nel_equi=length_tobe_covered_equi/len_el_th;
needed_equi=floor(nel_equi);

new_ratio=length_tobe_covered_equi/needed_equi;

x_mid=zeros(needed_equi+1,1);
x_mid(1)=xel(x_cheb_reliable);
for ii=1:needed_equi
    x_mid(ii+1)=x_mid(ii)-new_ratio;
end
nodes=[xel(1:x_cheb_reliable);x_mid;-flipud(x_mid);xel(end-x_cheb_reliable+1:end)];

nodes_r=unique(nodes)';

number_of_nodes=length(nodes)-3;
number_of_elements=number_of_nodes-1;

nodes_final=zeros(number_of_nodes,1);

rrr=0;
count=0;
for i=1:length(nodes_r)
    if abs(nodes_r(i))>1e-6
        count=count+1;
        nodes_final(count)=nodes_r(i);
    else
        if nodes_r(i)==0 && rrr==0
            count=count+1;
            rrr=1;
            nodes_final(count)=0;
        end
    end
end

nodes_final
number_of_nodes
number_of_elements

gll_nodes=lglnodes(N)/2+0.5;

delta_yw_plus=abs((nodes_final(1)-nodes_final(2)))*Re_tau*gll_nodes(end-1)
delta_yc_plus=(nodes(length(nodes)/2-1))*Re_tau*(gll_nodes(end/2)-gll_nodes(end/2+1))
y1_plus=0
y2_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-1)
y3_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-2)
y4_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-3)
y5_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-4)
y6_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-5)
y7_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-6)
y8_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-7)
% y9_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-8)
% y10_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-9)
% y11_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-10)
% y12_plus=abs(nodes_final(1)-nodes_final(2))*Re_tau*gll_nodes(end-11)

nel_z=nelz
delta_z_avg=Lz/((N+1)*nel_z)*Re_tau
delta_z_max=delta_z_avg*(gll_nodes(end/2)-gll_nodes(end/2+1))/(1/(N+1))
delta_z_min=delta_z_avg*gll_nodes(end-1)/(1/(N+1))

for i=1:number_of_elements
    dy_final(i)=nodes_final(i+1)-nodes_final(i);
end

for i=1:16
    dy_final2(i)=mesh_180_prev(i+1)-mesh_180_prev(i);
end

figure(1);
plot(nodes_final,ones(number_of_nodes),'ob'); hold on;
plot(nodes_final,[0 dy_final],'o-r'); hold on;

plot(mesh_180_prev,ones(17),'og'); hold on;
plot(mesh_180_prev,[0 dy_final2],'o-g'); hold on;

y_last=0;
y_last2=0;
for i=1:number_of_elements/2
    for j=1:N+1
        mesh_final((i-1)*(N+1)+j,1)=gll_nodes(end-j+1)*dy_final(i)*Re_tau+y_last;
        mesh_final2((i-1)*(N+1)+j,1)=gll_nodes(end-j+1)*dy_final2(i)*Re_tau+y_last2;
        if j==N+1 
            y_last=gll_nodes(end-j+1)*dy_final(i)*Re_tau+y_last;
            y_last2=gll_nodes(end-j+1)*dy_final2(i)*Re_tau+y_last2;
        end
    end
end

figure(2);
plot(mesh_final,'ob'); hold on;
plot(moser(1).yp,'or'); hold on;
plot(mesh_final2,'og'); hold on;

legend('Script from 360','Moser','Chebyshev (current)','Location','NorthWest');