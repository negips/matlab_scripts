%--------------------------------------------------------------------------
%generate GLL points
%--------------------------------------------------------------------------
%import gll points
% [b,hdr,tag,N1,nel,Lvar,sts] = readnek_mpi(2,'le','./small/meshsmall0.f00001','xup');
% N=N1-1;
% [nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(['small/','meshsmall0.f00001']);
[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(['saab_wing2d.IC']);
N1=lr1(1);
N=N1-1;
nel=length(elmap);
%Loop over all elements
% plotgll = zeros(length(ELsmall)*(N+1)^2,2);
plotgll = zeros(nel*(N+1)^2,2);

%for i = 1:length(ELsmall)
% for i = 1:nel
%     ELsmall(i).GLL(:,1) = squeeze(reshape(b{i}(:,:,3),N1^2,1,1));
%     ELsmall(i).GLL(:,2) = squeeze(reshape(b{i}(:,:,4),N1^2,1,1));
%     plotgll((i-1)*(N+1)^2+1:i*(N+1)^2,:) = ELsmall(i).GLL;
% end
for i = 1:nel
    EL2D(i).GLL(:,1) = squeeze(nekdata(i,:,1));
    EL2D(i).GLL(:,2) = squeeze(nekdata(i,:,2));
    plotgll((i-1)*(N+1)^2+1:i*(N+1)^2,:) = EL2D(i).GLL;
end
plot(plotgll(:,1),plotgll(:,2),'*','Markersize',2)
axis equal
%clear plotgll
%hold on
%plot(xpro,ypro,'k.-')
