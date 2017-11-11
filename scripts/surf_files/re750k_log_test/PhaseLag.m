function lsq = PhaseLag(par,time,cl,omega,alpha0,dalpha,xfoilalpha,xfoilcl)

  phi=par(1)*pi/180;          % phase lag quasi-steady
  intg_const=par(2);          % integration
%  ofst = par(3);
  theta = par(3);             % phase lag added-mass
  ofst = par(4);    

  phase_lag = -pi/2;          % initial phase

  alpha = alpha0 + dalpha*sin(omega*(time-6.0) + phase_lag);

  pitch = dalpha*sin(omega*(time-6.0) + phi + phase_lag);
  alpha_lagg = alpha0 + pitch;

  added_mass = abs(intg_const)*cos(omega*(time-6.0) + theta + phase_lag);

  cl_lagg = interp1(xfoilalpha,xfoilcl,alpha_lagg,'linear');

  cl_pred = ofst + added_mass + cl_lagg;

%  dbstop in PhaseLag at 16  
  lsq=sum((cl-cl_pred).^2);

end
