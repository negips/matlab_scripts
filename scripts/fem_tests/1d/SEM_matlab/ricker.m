function f = ricker(t,f0,t0)
arg = pi*f0*(t-t0);
arg = arg.*arg;
f = (2*arg-1).*exp(-arg);
