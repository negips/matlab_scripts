function [cz_pred] = rebuild_unsteady_force(time,cz,cz_model,alpha_model,k,uoo,alpha0,dalpha,theta,phi)

  c=0.5;
  b=c/2;
  omega = 2*k*uoo/c;

  t_C = besselh(1,2,k)./(besselh(1,2,k) + i*besselh(0,2,k));
  t_F=real(t_C);
  t_G=imag(t_C);
  omega_k = k*uoo/b;
  a=(0.5-0.15)/2/b;
  r2 = sqrt(1 + (a*k)^2);
  added_m_eff_amp = -pi*dalpha*k.*r2;
  added_m_phi = asin(1/r2);

  alpha = alpha0 + dalpha*sin(omega*time + theta + phi);
  dalphadt = dalpha*k*cos(omega*time + theta + phi);
  a=(0.35-0.5)/2/b;
  lagterm = dalphadt*(1/2-a);
  alpha_eff = alpha + lagterm;
  Ck = besselh(1,2,k)./(besselh(1,2,k) + i*besselh(0,2,k));
  alpha_eff_wake = real(alpha_eff*Ck);
  %r = sqrt(1 + 0.5-a);
  %ptch_amp_eff = ptch_amp*r*180/pi;

  cz_lagg = interp1(alpha_model,cz_model,alpha_eff*180/pi,'linear');

  added_mass = -added_m_eff_amp*sin(omega*time + theta + phi);

  cz_pred = added_mass + cz_lagg;

end
