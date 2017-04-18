% Testing 3d surface output

clear
clc
close all

ifsurf = 1;

Re=100E+3;
nu=1/Re;


for i=1:30
    fname = sprintf('%s%5.5d', 'vtksaab_wing0.f',i);
    
    [datastruc,data,lr1,nelt,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek_surf(fname, ifsurf);
    
%    figure(i)
%    scatter3(datastruc(1).data(:),datastruc(2).data(:),datastruc(3).data(:),[],datastruc(5).data(:), '.'); colorbar
%    view([-5 60])
%    legend(num2str(time))

    top=find(datastruc(5).data(:)>0);
    bot=find(datastruc(5).data(:)<0);
    
    X = datastruc(1).data(:);
    Y = datastruc(2).data(:);
    Z = datastruc(3).data(:);
    SNX = datastruc(4).data(:);
    SNY = datastruc(5).data(:);
    ut = nu*datastruc(6).data(:);
    p_ut = datastruc(7).data(:);
    
    XT=X(top);
    YT=Y(top);
    ZT=Z(top);
    p_utT=p_ut(top);
    utT = ut(top);
    
    %% remove some points
    ind = find(XT>0.25);
    tmpx=XT(ind);
    tmpy=YT(ind);
    tmpz=ZT(ind);
    tmput=utT(ind);
    tmpput=p_utT(ind);
    
    figure(1)
    scatter3(tmpx,tmpy,tmpz,[],tmput, '.'); 
    caxis([-0.01 0.02])   
    colorbar
    view([-5 60])
    title(num2str(time))

    pause(0.1)  

end%
