function [ x_unq indicies ind_unq n_unq] = real_unique( x, tol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
  tol=1e-12;
end  

tol=abs(tol);

[x2 isort]= sort(x);
x_unq = x2(1);
unq_c(1) = 1;
ind_unq(1) = 1;
n_unq(1) = 1;
indicies{1} = [1];
% disp(['Array size: ', num2str(length(x2))])
for i=2:length(x2);
   if (abs(x2(i)-x_unq(unq_c)>tol))
       indicies = [indicies [i]];
       unq_c=unq_c+1;
       x_unq = [x_unq x2(i)];
       ind_unq = [ind_unq i];
       n_unq = [n_unq 1];
   else
       n_unq(unq_c) = n_unq(unq_c)+1;
       indicies{unq_c} = [indicies{unq_c} i];
   end
   
end


end

