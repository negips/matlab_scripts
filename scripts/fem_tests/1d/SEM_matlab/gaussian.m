function f = gaussian(t,f0,t0)
if isempty(t0), t0 = 0.45/f0; end
f = exp(-( 2*pi*f0*(t-t0) ).^2) ;
