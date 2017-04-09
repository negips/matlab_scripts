function [smoothval]=smoothstep(xin,xmin,xmax)

  ind=find(xin<xmin);
  xin(ind) = xmin;
  ind=find(xin>xmax);
  xin(ind) = xmax;

  xfrac=(xin-xmin)/(xmax-xmin)*4 - 2;
  smoothval=(tanh(xfrac) - tanh(-2))./(tanh(2)-tanh(-2));             % taking tanh range from -2 to 2.
