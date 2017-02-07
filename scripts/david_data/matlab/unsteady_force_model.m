function lsq = unsteady_force_model(par,time,cz,cz_model,alpha_model,k,uoo,alpha0,dalpha,theta)

  phi=par(1);
  intg_const=par(2);
%  phi2=par(3);
  c=0.5;
  omega = 2*k*uoo/c;

  pitch = dalpha*sin(omega*time+theta);
  alpha_pred=alpha0+pitch;
  alpha_lagg=alpha0+dalpha*sin(omega*time+theta+phi);

  inst_OMEGA=omega*dalpha*cos(omega*time+theta);

  p_motion = intg_const*(inst_OMEGA);

  cz_lagg = interp1(alpha_model,cz_model,alpha_lagg*180/pi,'linear');

  cz_pred = p_motion + cz_lagg;
  
  lsq=sum((cz-cz_pred).^2);

end
