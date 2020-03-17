function [xu,yu,alphau,xl,yl,alphal]=naca_prof(c,x,t,m,p)

%Symmetric NACA profile
yt=5*t*c*(0.2969*sqrt(x/c)+(-0.1260)*x/c+(-0.3516)*(x/c).^2+0.2843*(x/c).^3+(-0.1036)*(x/c).^4);

yu0=yt;
yl0=-yt;

%Cambered NACA profile
ip=find(x<=p*c,1,'last');
yc(1:ip)=m*x(1:ip)/p^2.*(2*p-x(1:ip)/c);
yc(ip+1:length(x))=m*(c-x(ip+1:end))/(1-p)^2.*(1+x(ip+1:end)/c-2*p);

dycdx(1:ip)=2*m/p^2*(p-x(1:ip)/c);
dycdx(ip+1:length(x))=2*m/(1-p)^2*(p-x(ip+1:end)/c);

theta=atan(dycdx);

xu=x-yt.*sin(theta);
xl=x+yt.*sin(theta);
yu=yc+yt.*cos(theta);
yl=yc-yt.*cos(theta);

dytdx=5*t*(0.14845./(sqrt(x/c))-0.126-0.7032*(x/c)+0.8529*(x/c).^2-0.4144*(x/c).^3);

dyudx=dycdx+dytdx.*cos(theta);
dyldx=dycdx-dytdx.*cos(theta);

alphau=atan(dyudx);
alphal=atan(dyldx);

end

