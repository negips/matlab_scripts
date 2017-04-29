% Testing 3d surface output

clear
clc
close all

fname = 'vtksaab_wing0.f00001';
ifsurf = 1;

[datastruc,data,lr1,nelt,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek_surf(fname, ifsurf);

% figure
% scatter3(datastruc(1).data(:),datastruc(2).data(:),datastruc(3).data(:),[],datastruc(7).data(:), '.'); colorbar
% view([-5 60])
% legend(num2str(time))

top=find(datastruc(5).data(:)>0);
bot=find(datastruc(5).data(:)<0);

X = datastruc(1).data(:);
Y = datastruc(2).data(:);
Z = datastruc(3).data(:);
SNX = datastruc(4).data(:);
SNY = datastruc(5).data(:);
ut = datastruc(6).data(:);
p_ut = datastruc(7).data(:);

XT=X(top);
YT=Y(top);
ZT=Z(top);
p_utT=p_ut(top);
utT = ut(top);

%% remove some points
ind = find(XT>0.15);
tmpx=XT(ind);
tmpy=YT(ind);
tmpz=ZT(ind);
tmput=utT(ind);
tmpput=p_utT(ind);

%figure
%scatter3(tmpx,tmpy,tmpz,[],tmput, '.'); colorbar
%view([-5 60])
%legend(num2str(time))

%[zsort ind] = sort(ZT);
[xsort ind] = sort(XT);
ysort = YT(ind);
p_sep_sort = p_utT(ind);
zsort = ZT(ind);
var = utT(ind);

% unq_z = real_unique(ZT,1e-12);
[unq_x indicies ind_unq n_unq] = real_unique( XT, 1e-12 );

for i=1:length(unq_x)
  z_avg(i) = mean(p_sep_sort(indicies{i}));
end
figure
plot(unq_x,z_avg);

ind1 = min(n_unq);
ind2 = max(n_unq);









