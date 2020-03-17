%     Plot impulse response for inhomogenous Ginzberg Landau equation

clear
clc
close all

%f1='GZ_Uo4_muo4_mut_o0_mu1_-0.04_mut_conv0_mut_abs0_xofst5.mat';
%f2='GZ_Uo4_muo4_mut_o0_mu1_-0.04_mut_conv0_mut_abs0_xofst20.mat';
%f3='GZ_Uo4_muo4_mut_o0_mu1_-0.04_mut_conv0_mut_abs0_xofst40.mat';

f1='GZ_Uo4_muo4_mut_o0_mu1_-0.02_mut_conv0_mu_abs0_mut_abs0_xofst10.mat';
f2='GZ_Uo4_muo4_mut_o0_mu1_-0.02_mut_conv0_mu_abs0_mut_abs0_xofst50.mat';
f3='GZ_Uo4_muo4_mut_o0_mu1_-0.02_mut_conv0_mu_abs0_mut_abs0_xofst100.mat';

i1=load(f1);
i2=load(f2);
i3=load(f3);

tol=1e-1;

[r c] = size(i1.u_save);
i1.u_save2=i1.u_save;
for ii=1:c
  ind=find(i1.u_save(:,ii)<tol);
  i1.u_save2(ind,ii) = nan;
end

[r c] = size(i2.u_save);
i2.u_save2=i2.u_save;
for ii=1:c
  ind=find(i2.u_save(:,ii)<tol);
  i2.u_save2(ind,ii) = nan;
end

[r c] = size(i3.u_save);
i3.u_save2=i3.u_save;
for ii=1:c
  ind=find(i3.u_save(:,ii)<tol);
  i3.u_save2(ind,ii) = nan;
end

figure(1)
colormap('jet')
s1=surf(i1.xgll(:),i1.t_save/i1.Tosc,transpose(i1.u_save2),'EdgeColor','none');
hold on
s2=surf(i2.xgll(:),i2.t_save/i2.Tosc+0.5,transpose(i2.u_save2),'EdgeColor','none');
s3=surf(i3.xgll(:),i3.t_save/i3.Tosc+0.895,transpose(i3.u_save2),'EdgeColor','none');
colorbar
view(2)
ylim([0 1.5])
ylabel('time')
xlabel('x')

alpha(s1,0.2)
alpha(s2,0.4)
alpha(s3,0.6)

%SaveFig(gcf,'impulse_response_decreasing.png','plots/',1)

%% increasingly unstable

%f4='GZ_Uo4_muo0_mut_o0_mu1_0.04_mut_conv0_mu_abs0_mut_abs0_xofst20.mat';
%f5='GZ_Uo4_muo0_mut_o0_mu1_0.04_mut_conv0_mu_abs0_mut_abs0_xofst40.mat';
%f6='GZ_Uo4_muo0_mut_o0_mu1_0.04_mut_conv0_mu_abs0_mut_abs0_xofst60.mat';

f4='GZ_Uo4_muo0_mut_o0_mu1_0.02_mut_conv0_mu_abs0_mut_abs0_xofst10.mat';
f5='GZ_Uo4_muo0_mut_o0_mu1_0.02_mut_conv0_mu_abs0_mut_abs0_xofst50.mat';
f6='GZ_Uo4_muo0_mut_o0_mu1_0.02_mut_conv0_mu_abs0_mut_abs0_xofst100.mat';

i4=load(f4);
i5=load(f5);
i6=load(f6);

tol=1e-1;

[r c] = size(i4.u_save);
i4.u_save2=i4.u_save;
for ii=1:c
  ind=find(i4.u_save(:,ii)<tol);
  i4.u_save2(ind,ii) = nan;
end

[r c] = size(i5.u_save);
i5.u_save2=i5.u_save;
for ii=1:c
  ind=find(i5.u_save(:,ii)<tol);
  i5.u_save2(ind,ii) = nan;
end

[r c] = size(i6.u_save);
i6.u_save2=i6.u_save;
for ii=1:c
  ind=find(i6.u_save(:,ii)<tol);
  i6.u_save2(ind,ii) = nan;
end

figure(2)
colormap('jet')
s4=surf(i4.xgll(:),i4.t_save/i4.Tosc,transpose(i4.u_save2),'EdgeColor','none');
hold on
s5=surf(i5.xgll(:),i5.t_save/i5.Tosc+0.025,transpose(i5.u_save2),'EdgeColor','none');
s6=surf(i6.xgll(:),i6.t_save/i6.Tosc+0.075,transpose(i6.u_save2),'EdgeColor','none');
colorbar
view(2)
ylim([0 0.6])
ylabel('time')
xlabel('x')

alpha(s4,0.2)
alpha(s5,0.4)
alpha(s6,0.6)

%SaveFig(gcf,'impulse_response_increasing.png','plots/',1)


%% Highly unstable blob

%f7='GZ_Uo4_muo2.5_mut_o0_mu1_0_mut_conv0_mut_abs0_xofst20.mat';
%f8='GZ_Uo4_muo2.5_mut_o0_mu1_0_mut_conv0_mut_abs0_xofst55.mat';

f7= 'GZ_Uo4_muo1.5_mut_o0_mu1_0_mut_conv0_mu_abs2_mut_abs0_xofst20.mat';
f8= 'GZ_Uo4_muo1.5_mut_o0_mu1_0_mut_conv0_mu_abs2_mut_abs0_xofst50.mat';
f9= 'GZ_Uo4_muo1.5_mut_o0_mu1_0_mut_conv0_mu_abs2_mut_abs0_xofst100.mat';
f10= 'GZ_Uo4_muo1.5_mut_o0_mu1_0_mut_conv0_mu_abs2_mut_abs0_xofst5.mat';

i7=load(f7);
i8=load(f8);
i9=load(f9);
i10=load(f10);

%close(3)
%close(2)

tol=1e-1;

[r c] = size(i7.u_save);
i7.u_save2=i7.u_save;
for ii=1:c
  ind=find(i7.u_save(:,ii)<tol);
  i7.u_save2(ind,ii) = nan;
end

[r c] = size(i8.u_save);
i8.u_save2=i8.u_save;
for ii=1:c
  ind=find(i8.u_save(:,ii)<tol);
  i8.u_save2(ind,ii) = nan;
end

[r c] = size(i9.u_save);
i9.u_save2=i9.u_save;
for ii=1:c
  ind=find(i9.u_save(:,ii)<tol);
  i9.u_save2(ind,ii) = nan;
end

[r c] = size(i10.u_save);
i10.u_save2=i10.u_save;
for ii=1:c
  ind=find(i10.u_save(:,ii)<tol);
  i10.u_save2(ind,ii) = nan;
end

figure(3)
colormap('jet')
s7=surf(i7.xgll(:),i7.t_save/i7.Tosc,transpose(i7.u_save2),'EdgeColor','none');
hold on
s8=surf(i8.xgll(:),i8.t_save/i8.Tosc+0.068,transpose(i8.u_save2),'EdgeColor','none');
s9=surf(i9.xgll(:),i9.t_save/i9.Tosc+0.15,transpose(i9.u_save2),'EdgeColor','none');
%s10=surf(i10.xgll(:),i10.t_save/i10.Tosc,transpose(i10.u_save2),'EdgeColor','none');

colorbar
view(2)
ylabel('time')
xlabel('x')
ylim([0 0.35])

alpha(s7,0.2)
alpha(s8,0.4)
alpha(s9,0.6)
%alpha(s10,0.1)

%SaveFig(gcf,'impulse_response_separation_bubble.png','plots/',1)






