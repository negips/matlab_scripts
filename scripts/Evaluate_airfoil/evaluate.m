function [ result ] = evaluate( x,z,param,side )
% [ result ] = evaluate( x,z )
% This function computes pressure coefficient around an airfoil based on a
% a panel method, calculates the BL profiles, and perform stability
% analysis based on a database method.
%
% Input:
% x,y : normalized coordinates of airfoil in a normal to leading edge system
% param: structure including flow parameters 
%
% result.xn  : coordinates of x-stations
% result.nts : Envelope of N-factors for TS-waves
% result.ncf : Envelope of N-factors for CF-waves
% result.cp  : pressure coeeficient cp
% result.u   : inviscid streamwise velocity on surface of airfoil

command0=(['unset GFORTRAN_STDIN_UNIT; unset GFORTRAN_STDOUT_UNIT;']);
%
Re   = param.Re;
Mach = param.Mach;
Tref = param.Tref;
sweep = param.sweep;
alpha = param.alpha;
chord = param.chord;
fmax  = param.fmax;
%
[u,cp]=mypanel(x,z,alpha);
%
%
%transform to 3D
%
cp=cp*cosd(sweep)^2;
u=sqrt(1-cp-sind(sweep)^2);
[cpmax, istag]=max(cp);
u(istag)=0;
x=x*chord/cosd(sweep);
z=z*chord;

result.cp=cp;
result.u=u;

dd=[x z u x*0 x*0 (1:size(x,1))']';

if (side==1)
    icase=1001;
    dd=dd(:,istag:-1:1);
else
    icase=1002;
    dd=dd(:,istag:1:end);
end

%   generate geo.dat
filename=(['geo_u.dat_',num2str(icase)]);
fid = fopen(filename,'w');
fprintf(fid,'# chord   Re    Mach    Tinf    Sweepang\n');
fprintf(fid,'%13.6g %13.6g %13.6g %13.6g %13.6g\n',1,Re,Mach,Tref,sweep);
fprintf(fid,'# x/c      z/c     ue      cq      Tw      edge_index\n');
fprintf(fid,'%15.8g %15.8g %15.8g %15.8g %15.8g %i\n',dd );
fclose(fid);

%   generate bl.in
genblin( icase,length(x) )

%   run the BL code
command=([command0,'echo -e "',num2str(icase),'\nn\n" | ~/bin/bl3D > bl.out']);
system(command);
command=[command0,'rm geo_u.dat* geo.upd.dat* bl.in* stab.dat']; 
system(command);

%   run  database method
command=[command0,'~/bin/parabola ',num2str(fmax)];
system(command);
command=[command0,'rm histo.dat meanfield parabola.out']; 
system(command);

%   read envelope file
dd=importdata('envelop.dat',' ',3);
result.xn=dd.data(:,1);
result.ncf=dd.data(:,6);
result.nts=dd.data(:,8);

end

