function [ ] = genblin( icase,nx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
filename=(['bl.in_',num2str(icase)]);
fid = fopen(filename,'w');
text{1}='200                Number of points in normal direction';
text{2}='20                etamax';
text{3}='1                mutask (1:sutherland, 2:polynomial)';
text{4}='1                kptask (1:keye,       2:polynomial)';
text{5}='1                prtask (1:constant,   2:variable Prandtl number';
text{6}='1                cptask (1:constant,   2:variable cp)';
text{7}='0.720                Prandtl number: (0 to be computed)';
text{8}='1.400 287.03         Gamma and gas constant';
text{9}=([num2str(nx),'                Index of last station']);
text{10}=(['-1 2 ',num2str(nx),'   1      ftype, ifirst, ilast, istep']);
text{11}='0                Number of new points in each interval';
text{12}='1                Order of interpolation';
text{13}='0                isim';
text{14}='0                bcmtask: 0 mw read from geo, 1 pbox uses pbox.in';
text{15}='0 1 0 1 0 1          Output to files on/off';
text{16}='0                Curvature effects on/off (should always be off for now)';
text{17}='1000                index of starting position for turbulent computations';

for k=1:length(text)
    fprintf(fid,'%s\n',text{k});
end
fclose(fid);

end

