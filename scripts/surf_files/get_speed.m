% Get convective speed

clear
clc
close all


load('tr750k_n9_2.mat');

ind1=tr_time>48.44;
ind2=tr_time<49.34;
ind3=find(ind1.*ind2);

t=tr_time(ind3);
y=trx_uv(ind3);
c1 = polyfit(t,y,1)
y1 = polyval(c1,t);

ind4=tr_time>51.12;
ind5=tr_time<52.68;
ind6=find(ind4.*ind5);

t=tr_time(ind6);
y=trx_uv(ind6);
c2 = polyfit(t,y,1)
y2 = polyval(c2,t);


% Second loop
ind7=tr_time>56.26;
ind8=tr_time<57.33;
ind9=find(ind7.*ind8);

t=tr_time(ind9);
y=trx_uv(ind9);
c3 = polyfit(t,y,1)
y3 = polyval(c3,t);

ind10=tr_time>59.28;
ind11=tr_time<60.42;
ind12=find(ind10.*ind11);

t=tr_time(ind12);
y=trx_uv(ind12);
c4 = polyfit(t,y,1)
y4 = polyval(c4,t);


figure(1)
plot(tr_time,trx_uv); hold on

indt = [ind3(1) ind3(end) ind6(1) ind6(end) ind9(1) ind9(end) ind12(1) ind12(end)];
plot(tr_time(indt), trx_uv(indt), 's ', 'MarkerSize', 8, 'LineWidth', 2)


figure(2)
plot(tr_time(ind3),trx_uv(ind3));hold on
plot(tr_time(ind3),y1);hold on

plot(tr_time(ind6),trx_uv(ind6));
plot(tr_time(ind6),y2);hold on

plot(tr_time(ind9),trx_uv(ind9));hold on
plot(tr_time(ind9),y3);hold on

plot(tr_time(ind12),trx_uv(ind12));
plot(tr_time(ind12),y4);hold on



