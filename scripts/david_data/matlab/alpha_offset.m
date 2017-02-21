function lsq = alpha_offset(par,alpha,cz,alpha_model,cz_model)

  alpha_ofst=par(1);
  alpha_new = alpha+alpha_ofst;

  cz_model = interp1(alpha_model,cz_model,alpha_new*180/pi,'linear');
  
  lsq=sum((cz-cz_model).^2);

end
