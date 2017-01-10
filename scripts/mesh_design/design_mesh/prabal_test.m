% Testing _cheb/uniform mesh distribution
Nel = 2;
crit = 0;
lencrit = 0;

while crit==0

[xel, DM] = chebdif(Nel, 1);

%Distance between Chebychev nodes in element [-1,1]
d_line=abs(diff(xel));         
d_min_line=min(d_line)/2;
[d_max_line maxind]=max(d_line);
d_max_line = d_max_line/2;

cheb_scale = len_el_r/d_min_line;

scaled_d_max = cheb_scale*d_max_line;
tlength = cheb_scale*sum(d_line(1:maxind))/2;
max_arr=[max_arr scaled_d_max];

if len_el_th<scaled_d_max
     crit = 1;
     break;
end

if tlength > 1
     crit = 1;
     lencrit = 1;
     break;
end

Nel=Nel+1

end
Nel=Nel-1

