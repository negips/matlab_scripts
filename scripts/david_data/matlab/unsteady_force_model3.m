function lsq = unsteady_force_model3(par,time,cz,cz_model,alpha_model,k,uoo,alpha0,dalpha,theta)

  phi=par(1);
  intg_const=par(2);
  toff=par(3);
  steady_shift=par(4);
  c=0.5;
  omega = 2*k*uoo/c;

  time=time-toff;
  pitch = dalpha*sin(omega*time+theta);
  alpha_pred=alpha0+pitch;
  alpha_lagg=alpha0+dalpha*sin(omega*time+theta+phi);

  inst_OMEGA=omega*dalpha*cos(omega*time+theta);

  p_motion = intg_const*(inst_OMEGA);

  cz_lagg = interp1(alpha_model+steady_shift,cz_model,alpha_lagg*180/pi,'linear');

  cz_pred = p_motion + cz_lagg;
%  dbstop in unsteady_force_model at 21 
  
  lsq=sum((cz-cz_pred).^2);

end
