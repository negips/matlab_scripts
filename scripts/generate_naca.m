clear
clc
close all

iaf.designation='0012';
% designation='0008';
iaf.n=1000;
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
if (genfile)
  fid = fopen(['naca' iaf.designation '.dat'],'w');
  for i=length(af.xL):-1:1
    fprintf(fid,'%11.10f\t%11.10f\t%11.10f\n',af.xL(i),af.zL(i),0.);
    allx = [allx af.xL(i)];
    ally = [ally af.zL(i)];
  end 
  for i=length(af.xU)-1:-1:2
    fprintf(fid,'%11.10f\t%11.10f\t%11.10f\n',af.xU(i),af.zU(i),0.);
    allx = [allx af.xU(i)];
    ally = [ally af.zU(i)];
  end
  fclose(fid);
end



