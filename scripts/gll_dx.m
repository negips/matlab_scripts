% Change in dx

clear
clc
close all

n1=5;
n2=50;

N=n1:n2;

count=0;
for i=N
  count=count+1;
  [xgll wts] = lglnodes(i);
  dx(count) = min(abs(diff(xgll)));
end

plot(N,dx)
