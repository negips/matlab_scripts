function lsq = unsteady_alpha2(par,time,alpha,k,uoo)

  theta=par(1);
  alpha0=par(2);
  dalpha=par(3);
  c=0.5;
  
  omega = 2*k*uoo/c;
  
  alpha_pred=alpha0+dalpha*sin(omega*time+theta);
  
  lsq=sum((alpha-alpha_pred).^2);

end
