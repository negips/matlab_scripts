%     Plot impulse response for inhomogenous Ginzberg Landau equation

clear
clc
close all

f5='GZ_Uo4_muo4_mut_o0_mu1_-0.04_mut_conv0_mut_abs0_xofst5.mat';
f20='GZ_Uo4_muo4_mut_o0_mu1_-0.04_mut_conv0_mut_abs0_xofst20.mat';
f40='GZ_Uo4_muo4_mut_o0_mu1_-0.04_mut_conv0_mut_abs0_xofst40.mat';

i5=load(f5);
i20=load(f20);
i40=load(f40);

tol=1e-1;

[r c] = size(i5.u_save);
i5.u_save2=i5.u_save;
for ii=1:c
  ind=find(i5.u_save(:,ii)<tol);
  i5.u_save2(ind,ii) = nan;
end

[r c] = size(i20.u_save);
i20.u_save2=i20.u_save;
for ii=1:c
  ind=find(i20.u_save(:,ii)<tol);
  i20.u_save2(ind,ii) = nan;
end

[r c] = size(i40.u_save);
i40.u_save2=i40.u_save;
for ii=1:c
  ind=find(i40.u_save(:,ii)<tol);
  i40.u_save2(ind,ii) = nan;
end



figure(1)
colormap('jet')
s5=surf(i5.xgll(:),i5.t_save/i5.Tosc,transpose(i5.u_save2),'EdgeColor','none');
hold on
s20=surf(i20.xgll(:),i20.t_save/i20.Tosc+0.27,transpose(i20.u_save2),'EdgeColor','none');
s40=surf(i40.xgll(:),i40.t_save/i40.Tosc+0.5,transpose(i40.u_save2),'EdgeColor','none');
colorbar
view(2)

alpha(s5,0.2)
alpha(s20,0.4)
alpha(s40,0.6)

%SaveFig(gcf,'impulse_response.png','plots/',1)

%% increasingly unstable

h20='GZ_Uo4_muo0_mut_o0_mu1_0.04_mut_conv0_mu_abs0_mut_abs0_xofst20.mat';
h40='GZ_Uo4_muo0_mut_o0_mu1_0.04_mut_conv0_mu_abs0_mut_abs0_xofst40.mat';
h60='GZ_Uo4_muo0_mut_o0_mu1_0.04_mut_conv0_mu_abs0_mut_abs0_xofst60.mat';

ii20=load(h20);
ii40=load(h40);
ii60=load(h60);

tol=1e-1;

[r c] = size(ii20.u_save);
ii20.u_save2=ii20.u_save;
for ii=1:c
  ind=find(ii20.u_save(:,ii)<tol);
  ii20.u_save2(ind,ii) = nan;
end

[r c] = size(ii40.u_save);
ii40.u_save2=ii40.u_save;
for ii=1:c
  ind=find(ii40.u_save(:,ii)<tol);
  ii40.u_save2(ind,ii) = nan;
end

[r c] = size(ii60.u_save);
ii60.u_save2=ii60.u_save;
for ii=1:c
  ind=find(ii60.u_save(:,ii)<tol);
  ii60.u_save2(ind,ii) = nan;
end



figure(4)
colormap('jet')
ss20=surf(ii20.xgll(:),ii20.t_save/ii20.Tosc,transpose(ii20.u_save2),'EdgeColor','none');
hold on
ss40=surf(ii40.xgll(:),ii40.t_save/ii40.Tosc+0.03,transpose(ii40.u_save2),'EdgeColor','none');
ss60=surf(ii60.xgll(:),ii60.t_save/ii60.Tosc+0.08,transpose(ii60.u_save2),'EdgeColor','none');
colorbar
view(2)

alpha(ss20,0.2)
alpha(ss40,0.4)
alpha(ss60,0.6)

%SaveFig(gcf,'impulse_response.png','plots/',1)




%% Highly unstable blob

g20='GZ_Uo4_muo2.5_mut_o0_mu1_0_mut_conv0_mut_abs0_xofst20.mat';
g55='GZ_Uo4_muo2.5_mut_o0_mu1_0_mut_conv0_mut_abs0_xofst55.mat';

j20=load(g20);
j55=load(g55);

%close(3)
%close(2)

tol=1e-1;

[r c] = size(j20.u_save);
j20.u_save2=j20.u_save;
for ii=1:c
  ind=find(j20.u_save(:,ii)<tol);
  j20.u_save2(ind,ii) = nan;
end

[r c] = size(j55.u_save);
j55.u_save2=j55.u_save;
for ii=1:c
  ind=find(j55.u_save(:,ii)<tol);
  j55.u_save2(ind,ii) = nan;
end


figure(3)
colormap('jet')
t20=surf(j20.xgll(:),j20.t_save/j20.Tosc,transpose(j20.u_save2),'EdgeColor','none');
hold on
t55=surf(j55.xgll(:),j55.t_save/j55.Tosc+0.55,transpose(j55.u_save2),'EdgeColor','none');
colorbar
view(2)

alpha(t20,0.2)
alpha(t55,0.4)

%SaveFig(gcf,'impulse_response2.png','plots/',1)






