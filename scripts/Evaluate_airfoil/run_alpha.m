clear all
close all
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
setenv('PATH', ['/usr/bin:/bin:/usr/sbin:/sbin:.:/Users/ardeshir/bin:/opt/local/bin'])
%%
%
Pref=1e5;
Tref=300;
Rgas=287.15;
Rhoref=Pref/(Tref*Rgas);
Visc=1.846e-5;
Uref=30;
Lref=1.0;
Mach=Uref/sqrt(Rgas*1.4*Tref);
Re=Uref*Rhoref*Lref/Visc;
%%
%
myfile='asu67.dat';
xy=importdata(myfile);
%xy=flipud(xy);
%%
param.Re = Re  ;
param.Mach = Mach;
param.Tref = Tref;
param.fmax = 3e-4/(2*pi*Visc/Rhoref/Uref^2); % just a guess
str='';
k=0;
%%
side=1;
for chord=1 %0.5:0.1:0.7
    for sweep=35%:10:45
        for alfa_3D=-2:1:0
            
            param.alpha=atand(tand(alfa_3D)/cosd(sweep));
            param.sweep = sweep;
            param.chord = chord;
            
            k=k+1;
            a(k)=alfa_3D;
            s(k)=sweep;
            c(k)=chord;
            
            result=evaluate(xy(:,1),xy(:,2),param,side);
            Ncf(k)=max(result.ncf);
            Nts(k)=max(result.nts);
            figure(1)
            plot(xy(:,1),result.cp)
            hold on
            figure(2)
            plot(result.xn,result.ncf)
            hold on
            str1=sprintf('C=%g, Sweep=%g, AoA=%g \n',chord,sweep,alfa_3D);
            str=strcat(str,str1);
            %legend(str)
        end
    end
end
%%
%
mylegend='';
figure(3);clf;
for k=1:length(a)
plot(a(k),Ncf(k),'k+-',a(k),Nts(k),'ro-')
mylegend=([mylegend, sprintf('Ncf, C=%g, Sweep=%g',c(k),s(k)),'']);
end
xlabel('AoA','FontSize',18)
ylabel('Max N','FontSize',18)
title('ASU(67)-0315 profile','FontSize',18)
legend(mylegend)
grid on
