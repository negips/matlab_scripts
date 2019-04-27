function Ufsc=similar_fsc(y,mfac)
% Ufsc=SIMILAR_FSC(YIN,MFAC)
% This function computes the selfsimilar Falkner-Skan-Cook boundary-layer profiles
%
% Input
%	y:  contains normal coordinates weher velocity values are required.
%	mfac:  Falkner-Skan-Cooke parameter U=C*x^mfac. 
%
% Output 
%	Ufsc.u, Ufsc.v, Ufsc.w: streamwise, normal and spanwise velocity 
%	Ufsc.uy, Ufsc.vy, Ufsc.wy:  normal derivatives of u, v, and w
%
% (c) Ardeshir Hanifi & David Tempelmann
%
%-------------------------------------------------------------------------------
% Start approximation for f'' and g' 
fbis= -0.0791*mfac^4+0.7414*mfac^3-1.4302*mfac^2+1.6813*mfac+0.3318;
fbis=fbis*2/(mfac+1);
gprim= -1.3429*mfac^4+2.6745*mfac^3-1.8293*mfac^2+0.7102*mfac+0.3244;
gprim=gprim*2/(mfac+1);

% Options for ODE solver
options = odeset('RelTol',1e-12,'AbsTol',[1e-8 1e-8 1e-8 1e-8 1e-8]);
tspan=[0 max(y)];
ic=[0 0 fbis 0 gprim];

% Integrate equations until convergend solution is found
[yt,f]=ode45(@(t,y) fsc(t,y,mfac),tspan,ic,options);
ue=f(end,2);
ue1=ue;
fbis1=fbis;
fbis=fbis+0.01;
tol=1e-8;
while abs(ue-1)>tol
    ic=[0 0 fbis 0 1];
    [yt,f]=ode45(@(t,y) fsc(t,y,mfac),tspan,ic,options);
    ue=f(end,2);
    dfbis=(1-ue)*(fbis-fbis1)/(ue-ue1);
    fbis1=fbis;
    ue1=ue;
    fbis=fbis+dfbis;
end

% Interpolate the velocity profiles on the input grid
c=1/sqrt((mfac+1)/2);
vv=((1.0-mfac)*f(:,2).*yt-(1.0+mfac)*f(:,1))*.5*c;
vvz=((1.0-mfac)*f(:,2)+(1.0-mfac)*f(:,3).*yt-(1.0+mfac)*f(:,2))*.5*c;
Ufsc.u=interp1(yt*c,f(:,2),y,'pchip');
Ufsc.w=interp1(yt*c,f(:,4)/f(end,4),y,'pchip');
Ufsc.v=interp1(yt*c,vv,y,'pchip');
Ufsc.uy=interp1(yt*c,f(:,3)/c,y,'pchip');
Ufsc.wy=interp1(yt*c,f(:,5)/f(end,4)/c,y,'pchip');
Ufsc.vy=interp1(yt*c,vvz/c,y,'pchip');





    
