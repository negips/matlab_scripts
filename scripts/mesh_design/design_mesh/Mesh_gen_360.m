

close all; clear all; clc
format long
Re_tau=360;     % Reynolds number based on the friction velocity
N=11;            % Polynomial order. (Number of points on the element -1)
Nel=31;         % Number of elements to be put on the line (based on the chebychev grid) -1
%19
[x]=lglnodes(N);  % x represens GLL Nodes



%%%#######################################################################
%%%################  Grid in the Axial Direction (Z)  ####################
%%%#######################################################################
Lz=25;          % Length of the pipe. Based on this length, and Re_tau and also the number
% GLL nodes on each element, we can judge how many stages we need to cope
% with the mesh criteria in the z direction (A.K.A.dMaxElPlus=10 )
dz=abs(x(2:end)-x(1:end-1));
dz_max=max(dz)/2;

dzMaxElPlus=10; %z_plus<=10
dz_max_el=dzMaxElPlus/Re_tau;
len_el_z=dz_max_el/dz_max;
nelz=ceil(Lz/len_el_z) % number of stages that is needed to maintain 
%the axia grid criteria (z_plus<=10). This will be applied in n2to3 as
%number os slabes (stages).


%%%#######################################################################
%%%################  Grid in the wall_normal Direction (h)  ###################
%%%#######################################################################

dr=abs(x/2-0.5);
dr_min=dr(4);  % to maintain four points below rPluas=1.

drWallElPlus=1.0; %r_plus=1 to be maintained as a key criteria
dr_min_el=drWallElPlus/Re_tau;
len_el_r=dr_min_el/dr_min

dthMaxElPlus=4.5; %r_plus<=5 Touch this value to get the right behavior!!
dth=abs(x(2:end)-x(1:end-1));
dth_max=max(dth)/2;
dth_max_el=dthMaxElPlus/Re_tau;
len_el_th=dth_max_el/dth_max

% Now play with Nel to maintain the condition for rPlus
[xel, DM] = chebdif(Nel, 1);

d_line=abs(xel(2:end)-xel(1:end-1)); % distance between Chebychev nodes in element [-1,1]
d_min_line=min(d_line)
d_max_line=max(d_line)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% poly order 11  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Let's say the first part is Chebychev and the second part is covered by
%%% equi-distance grid blocks. In that case:
x_cheb_reliable=5
length_tobe_covered_equi=xel(x_cheb_reliable)
nel_equi=length_tobe_covered_equi/len_el_th
needed_equi=floor(nel_equi)

new_ratio=length_tobe_covered_equi/needed_equi

x_mid=zeros(needed_equi+1,1);
x_mid(1)=xel(x_cheb_reliable);
for ii=1:needed_equi
x_mid(ii+1)=x_mid(ii)-new_ratio;
end
nodes=[xel(1:x_cheb_reliable);x_mid;-flipud(x_mid);xel(end-x_cheb_reliable+1:end)];
figure ; plot(nodes,ones(length(nodes)),'r*')

nodes_final=unique(nodes)'

number_of_nodes=length(nodes)-3
number_of_elements=number_of_nodes-1
delta_yw_plus=(nodes(1)-nodes(2))*360*0.027550363888559
delta_yc_plus=(nodes(length(nodes)/2-1))*360*0.136552932854928
y1_plus=0
y2_plus=(nodes(1)-nodes(2))*360*0.027550363888559
y3_plus=(nodes(1)-nodes(2))*360*0.090360339177997
y4_plus=(nodes(1)-nodes(2))*360*0.183561923484070
y5_plus=(nodes(1)-nodes(2))*360*0.300234529517326
y6_plus=(nodes(1)-nodes(2))*360*0.431723533572536
y7_plus=(nodes(1)-nodes(2))*360*0.568276466427464
