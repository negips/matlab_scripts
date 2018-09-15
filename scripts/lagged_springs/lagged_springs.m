% Solve a multiple spring mass system
% Such that the different springs acting with a phase lag.

clear
clc
close all

K = [-1.0  -0.02];        % Spring stiffness constants
tlag = [-0.0  -0.01];       % time lag for each spring.

lagstart = 10000;

D = [-0.01];      % Dissipation

nK = length(K);

ifp = find(tlag>0.);
if ~isempty(ifp)
   disp(['Time Lags must be non positive'])
   tlag
   return
end

nsteps = 1000000;
 

% Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0 1  0 0];
ex1   = [0 2 -1 0];
ex2   = [0 3 -3 1];

x0 = 0;
v0 = 0.01;

xlag = zeros(3,1);
vlag = zeros(3,1);

vlag(1) = v0;
xlag(1) = x0;

istep = 0;
dt = 0.001;
time = 0.;

iostep = 5000;
if (iostep>0)
  figure(1)
end  

xhis=zeros(nsteps+1,1);
x2his=zeros(nsteps+1,1);
vhis=zeros(nsteps+1,1);
this=zeros(nsteps+1,1);


while istep<=nsteps
  
  istep = istep+1;
  time = time+dt; 

  % set time integration operators
  if istep == 1
       bdfk = bdf1;
       extk = ex0;
  elseif istep == 2
       bdfk = bdf2;
       extk = ex1;
  else
       bdfk = bdf3;
       extk = ex2;
  end

  b = zeros(2,1);
  bx = 0; bdx = 0;
  bv = 0; bdv = 0;
  x2 = 0;
  for i=1:nK
%   Include the spring only if time is greater than the time lag
%   Ensure time lags are all negative.
%   Phase gains can be modelled with large phase lags
    if tlag(i)==0
      tmp = K(i)*extk(2:end)*xlag(:,i)*dt;
      bv = bv + tmp;
    elseif (time+tlag(i))/abs(tlag(i))>lagstart
      tmp = K(i)*interp1(this(1:istep),xhis(1:istep),time+tlag(i),'pchip');
      bv = bv+tmp;
      x2 = tmp/K(i);
    end

  end  
  bx = extk(2:end)*vlag*dt;

  bdv = bdfk(2:end)*vlag;
  bdx = bdfk(2:end)*xlag;

%  dissip = D*extk(2:end)*vlag*dt;

  vn = [bv - bdv]/(bdfk(1)-D*dt);
  xn = [bx - bdx]/bdfk(1);

  vlag(2:end) = vlag(1:end-1);
  xlag(2:end) = xlag(1:end-1);
  vlag(1) = vn;
  xlag(1) = xn;

  xhis(istep) = xn;
  x2his(istep) = x2;
  vhis(istep) = vn;
  this(istep) = time;

  if (iostep>0 && mod(istep,iostep)==0) || istep==nsteps
    figure(1)
%    plot(this(1:istep),xhis(1:istep))
    plot(xhis(1:istep),vhis(1:istep))  
    pause(0.01)
  end

end    


[pks,locs]=findpeaks(xhis);
pktimes=this(locs);
figure(2)
plot(diff(pktimes))
      

        



