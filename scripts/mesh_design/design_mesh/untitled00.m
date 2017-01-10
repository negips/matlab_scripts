

close all; clear all; clc
format long
Re_tau=180;     % Reynolds number based on the friction velocity
N=11;            % Polynomial order. (Number of points on the element -1)
Nel=18;         % Number of elements to be put on the line (based on the chebychev grid) -1

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

dthMaxElPlus=5; %r_plus<=5
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
%%%%%%%%%%%%%%%%% poly order 7  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_ratio=0.587785252292473/5;
for ii=1:10
x_mid(ii)=0.587785252292473-new_ratio*ii;
end
nodes=[xel(1:7);x_mid';xel(16:end)];
figure ; plot(nodes,ones(length(nodes)),'r*')

