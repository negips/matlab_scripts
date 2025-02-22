clear
clc
close all

iaf.designation='0012';
% designation='0008';
iaf.n=100000;
iaf.HalfCosineSpacing=1;
iaf.wantFile=0;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;

af = naca4gen(iaf);

% plot(af.x,af.z,'bo-')

plot(af.xU,af.zU,'bo-')
hold on
plot(af.xL,af.zL,'ro-')

%axis equal
genfile=1;
allx = [];
ally = [];
fmt = '%16.14f\t%16.14f\t%16.14f\n';

if (genfile)
  fid = fopen(['naca' iaf.designation '.dat'],'w');
  for i=length(af.xL):-1:1
    fprintf(fid,fmt,af.xL(i),af.zL(i),0.);
    allx = [allx af.xL(i)];
    ally = [ally af.zL(i)];
  end 
  for i=length(af.xU)-1:-1:2
    fprintf(fid,fmt,af.xU(i),af.zU(i),0.);
    allx = [allx af.xU(i)];
    ally = [ally af.zU(i)];
  end
  fclose(fid)
end

% fid=fopen('add_curve.rpl');
% fid2 = fopen('add_new_curve_naca.rpl','w');
% 
% tline = fgetl(fid);
% npts = length(allx); 
% while ischar(tline)
%   ind = findstr(tline, '{pnt_fine0'); 
%   if isempty(ind)
%     fprintf(fid2,'%s\n',tline);
%   else
%     tline2 = tline(1:ind);  
%     for i=0:npts-1
%        tline2 = [tline2 'pnt_vv' num2str(i) ' '];
%     end
%     tline2 = [tline2 ' pnt_vv0}'];
%     fprintf(fid2,'%s\n', tline2);
%   end
%   tline = fgetl(fid);
% end

fclose all; 



