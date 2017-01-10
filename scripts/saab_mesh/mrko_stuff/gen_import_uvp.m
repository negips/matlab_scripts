%--------------------------------------------------------------------------
%generate GLL points
%--------------------------------------------------------------------------
%% import gll points from baseflow file
% - singleprecision file
% [b,hdr,tag,N1,nel,Lvar,sts] = readnek_mpi(2,'le','./base/wing25D0.f00112','xupt');
% N=N1-1;
% - double precision file(latest)
% [b,hdr,tag,N1,nel,Lvar,sts] = readnek_mpi(2,'le',importfile,'xup');
% N=N1-1;
% - double precision
clear
clc

[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(['saab_wing2d.IC']);%
N1=lr1(1);
N=N1-1;

%% squeeze data in element structure
for i = 1:length(elmap)
    EL2D(i).GLL(:,1) = squeeze(nekdata(i,:,1));
    EL2D(i).GLL(:,2) = squeeze(nekdata(i,:,2));
    EL2D(i).VEL(:,1) = squeeze(nekdata(i,:,3));
    EL2D(i).VEL(:,2) = squeeze(nekdata(i,:,4));
    EL2D(i).P(:,1) = squeeze(nekdata(i,:,5));
%    ELbase(i).VEL(:,3) = squeeze(nekdata(i,:,6));

%    scatter3(EL2D(i).GLL(:,1),EL2D(i).GLL(:,2),EL2D(i).VEL(:,1));colorbar; 
end

%% reshape element structure data in input form
gll2Dbase = zeros(length(EL2D)*(N+1)^2,2);
for i = 1:length(EL2D)
    gll2Dbase((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL2D(i).GLL(1:(N+1)^2,1:2);
    vel2Dbase((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL2D(i).VEL(1:(N+1)^2,1:2);
    p2Dbase((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL2D(i).P(1:(N+1)^2,1);
end


%% readnek
%   This function reads binary data from the nek5000 file format
%  
%     [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(fname)
%  
%     INPUT
%     - fname:  name of the file 
%  
%     OUTPUT
%     - data:   nek5000 data ordered as (iel,inode,[x|y|(z)|u|v|(w)|p|T|s_i])
%     - lr1:    element-size vector (lx1,ly1,lz1)
%     - elmap:  reading/writing map of the elements in the file
%     - time:   simulation time
%     - istep:  simulation step
%     - fields: fields saved in the file
%     - emode:  endian mode 'le' = little-endian, 'be' = big-endian
%     - wdsz:   single (4) or double (8) precision
%     - etag:   tag for endian indentification
%     - header: header of the file (string)
%     - status: status (< 0 something went wrong)
