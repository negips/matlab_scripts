%----------------------------------------------------------------------
%extend cut baseflow to 3D
%----------------------------------------------------------------------
gen_import_uvp

disp('Extrude')
% addpath('/scratch/mrko/work/matlab-library/nek')
%% Get GLL points                                                 
% - read cut baseflow file
% [nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(['small3D/',casename3D,'.IC']);
[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(['saab_wing3d0.f00001']);

nel = length(elmap);

% - save in EL structure  3D
for iel = 1:nel
    EL3D(iel).GLL(:,1) = squeeze(nekdata(iel,:,1));
    EL3D(iel).GLL(:,2) = squeeze(nekdata(iel,:,2));
    EL3D(iel).GLL(:,3) = squeeze(nekdata(iel,:,3));
end

nelz = nel/length(EL2D);
disp(['Spanwise elements: ' num2str(nelz)]);

for ielz = 1:nelz
    for iel2D = 1:length(EL2D)
        for k=1:lr1(3)
            for j=1:1
               iel = iel2D + length(EL2D)*(ielz-1);

%     ELsmall3D(iel).VEL(:,1) = ELsmall(iel2D).VEL(:,1);
%     ELsmall3D(iel).VEL(:,2) = ELsmall(iel2D).VEL(:,2);
% %     ELsmall3D(iel).VEL(:,3) = ELsmall(iel2D).VEL(:,3);
%     ELsmall3D(iel).P = ELsmall(iel2D).P;
%     
    EL3D(iel).VEL((k-1)*lr1(3)^2+1:k*lr1(3)^2,1) = EL2D(iel2D).VEL((j-1)*lr1(3)^2+1:j*lr1(3)^2,1);
    EL3D(iel).VEL((k-1)*lr1(3)^2+1:k*lr1(3)^2,2) = EL2D(iel2D).VEL((j-1)*lr1(3)^2+1:j*lr1(3)^2,2);
    EL3D(iel).VEL((k-1)*lr1(3)^2+1:k*lr1(3)^2,3) = 0*EL2D(iel2D).VEL((j-1)*lr1(3)^2+1:j*lr1(3)^2,1);
    EL3D(iel).P((k-1)*lr1(3)^2+1:k*lr1(3)^2,1) = EL2D(iel2D).P((j-1)*lr1(3)^2+1:j*lr1(3)^2,1);
            end
        end
    end
end


%% Write 3D initial condition
% reshape data for writing
nekdata = zeros(length(EL3D),lr1(1)*lr1(2)*lr1(3),8);%nel
for iel = 1:length(EL3D)%nel
    nekdata(iel,:,1) = EL3D(iel).GLL(:,1);
    nekdata(iel,:,2) = EL3D(iel).GLL(:,2);
    nekdata(iel,:,3) = EL3D(iel).GLL(:,3);
    nekdata(iel,:,4) = EL3D(iel).VEL(:,1);
    nekdata(iel,:,5) = EL3D(iel).VEL(:,2);
    nekdata(iel,:,6) = EL3D(iel).VEL(:,3);
    nekdata(iel,:,7) = EL3D(iel).P(:,1);
end
fields  = 'XUP';%T';
%% write initial condition file
casename3D='saab_wing3d';
disp('Write 3D initial condition')
%status = writenek([casename3D,'.IC'],nekdata,lr1,elmap,0,0,fields,emode,8,etag);%8,etag);
fname=[casename3D, '.IC'];

%hdr = '#std 4        1                          0.00000000000E+00         0      0      1 XUP                                              ';
hdr = char(zeros(1,132));
wdsz=4;
time=0.;
istep=0;
fidd=0;
nf=1;
fields='XUP';
hdr2 = sprintf('#std %1i %2i %2i %2i %10i %10i %20.13E %9i %6i %6i %s\n',...
                 wdsz,lr1(1),lr1(2),lr1(3),nel,nel,time,istep,fidd,nf,fields);

l1 = length(hdr2);
hdr(1:l1) = hdr2;
display(hdr)

tag=6.54321;
dim=3;
endian='le';
fields='xup';

status = schreib_sp(nekdata,nel,hdr,tag,1,dim,endian,fname,fields);
disp('Done!')

% save in simulation folder
% system(sprintf('cp %s.rea ../base-torun/%s',simname3D,simname3D));
% system(sprintf('cp %s.IC  ./small3D/%s',casename3D));
% system(sprintf('cp %s.IC  ./small3D/%s','6mm6z'));
