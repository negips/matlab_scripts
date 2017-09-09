% Complex wavenumber dispersion analysis
% following Thompson and Pinsky (1994)
% + discrete boundary impedance analysis
% for the 1D wave equation solved by SEM

h = 1;		% element size
rho = 1;	% density
mu = 1;		% shear modulus
c = sqrt(mu/rho);	% shear wave velocity

Ps = [7:8];	% polynomial degrees to be analyzed
KIND = 'gll'; 	% 'pro'=prolate; 'gll'= Legendre SEM
ADD_LAYER = 1;	% add a viscous/PML/other element next to boundary:
		%	1 = Kelvin-Voigt homogeneous
		%	2 = Kelvin-Voigt tapered
		%	3 = frequency-shifted Perfectly Matched Layer
ETA = 0.3;	% Kelvin-Voigt viscosity / grid-critical-dt
PML_A = 10;	% PML damping amplitude
PML_N = 4;	% PML damping decay order
PML_RWC = 0.2;	% PML cut-off frequency / grid-critical-frequency

% figures: 1= dispersion curve
%          2= max frequency
%          3= impedance
FIGS = [2 0 1]; 

np = length(Ps);
dxm = 4*h./((Ps+1).^2-1);	% GLL min dx (approx)
wmax = 7/3 * c./ dxm; 		% GLL max frequency (approx)

%-- loop on polynomial degree
for p=1:np,

  P = Ps(p); % polynomial degree
  NGLL = P+1; % number of GLL nodes per element
  
  [xgll,wgll,H] = GetGLL(NGLL,KIND);
  dx_dxi = 0.5*h;
  dxm_num(p) = (xgll(2)-xgll(1))*dx_dxi; % GLL min dx
  
  % build the elementary mass matrix
  M = rho*wgll*dx_dxi;
  
  % build the elementary stiffness matrix
  W = mu*wgll/dx_dxi;
  K = H * ( repmat(W,1,NGLL).* H');
  K = (K+K')/2; % chop round-off errors by enforcing symmetry
  
  % boundary coefficient
  B = 1;
  
  % stuff for viscous layer
  x1 = 0.5*(1+xgll);
  switch ADD_LAYER
    case 2, 
      eta_taper = exp(-pi*x1.^2); 
    case 3, 
      PML_ALPHA = PML_A*c/h*( 1-(x1/h) ).^PML_N; 
      PML_WC = PML_RWC*wmax(p);
  end
  
  % frequencies (omega) to be analyzed:
  w = linspace(0.1, 2*wmax(p), 1000);
  %w = logspace(0, 2, 1000);
  
  k = zeros(length(w),1);	% wavenumber(omega)
  Z = ones(length(w),1);	% impedance(omega)

  ib = [1 NGLL]; 	% boundary nodes
  ii = [2:NGLL-1]; 	% interior nodes

  % stuff for fixing wrap-around problem
  rkpre = 0; shift = 0; factor = 1;

  % start loop on frequencies ...
  for iw = 1:length(w),

   %-- DISPERSION RELATION
  
    % elementary dynamic matrix :
    %S = -w(iw)^2 *diag(M) + K;
          iom = 1i*w(iw); % i*omega
          S = iom^2 *diag(M) + (1+ETA*dtc*iom)*K;
    % condensed on the boundary nodes:
    tmp = S(ii,ii)\S(ii,ib);
    g = S(ib,ib) - S(ib,ii)*tmp;
    % With the following notations:
    G1 = g(1,2);
    G2 = -(g(1,1)+g(2,2))/2;
    % the discrete problem (only boundary nodes) becomes:
    %   G1*d[n-1] -2*G2*d[n] +G1*d[n+1] =0
    % Assuming d[n] prop.to exp(i*k*h*n) the dispersion relation is obtained:
    kh = acos(G2/G1);
  
    % try to fix wrap-around
      khbis = factor*kh +shift;
      if real(khbis)<rkpre, 
        shift=shift+2*pi; 
        factor=-factor; 
        khbis = factor*kh +shift;
      end
      rkpre = real(khbis);
    % ... real part still needs unwrapping, see below
  
    % fix decay
      if imag(khbis)<0, khbis=conj(khbis); end
  
    k(iw) = khbis/h;
  
   %-- BOUNDARY IMPEDANCE
  
    % the discrete problem on the terminal element is
    %   -G2*d[1] +G1*d[2] = -B*T
    % hence
    %   -G2 + G1*exp(i*k*h) = -B*T
    % Combining with the dispersion relation, we obtain the
    % boundary node impedance (Z=T/v) :  
    %   Z = sin(k*h)/omega *G1/B
  
    %Z(iw) = sin(khbis)/w(iw) *G1/B;
    Z(iw) = sqrt(G1^2-G2^2)/w(iw) /B;
  
    % add a 1-element layer of some viscous material
    % and compute the new boundary impedance
    if ADD_LAYER
      iom = 1i*w(iw); % i*omega
      dtc = 2/wmax(p); % critical timestep (undamped)
      switch ADD_LAYER
        case 1, % Kelvin-Voigt uniform
          Sl = iom^2 *diag(M) + (1+ETA*dtc*iom)*K;
        case 2, % Kelvin-Voigt tapered
          Sl = iom^2 *diag(M) + K +ETA*dtc*iom*K*diag(eta_taper);
        case 3,  
          Gamma = (iom+PML_WC)./(iom+PML_WC+PML_ALPHA);
          Sl = iom^2 *diag(M./Gamma) + H*diag(Gamma.*W)*H' ;
          %Sl = iom^2 *diag(M) + diag(Gamma)*H*diag(Gamma.*W)*H' ;
      end
  
     % gl = condensed elementary matrix of the layer problem
      tmp = Sl(ii,ii)\Sl(ii,ib);
      gl = Sl(ib,ib) - Sl(ib,ii)*tmp;
     % coupling the layer to the infinite domain (index 2 for the interface):
     %   [gl]*(d1,d2) = (T1,-T2)
     %   T2 = i*omega*Z(infinite)*d2
     % leads to:
      Z(iw) = ( gl(1,1) -gl(1,2)*gl(2,1)/(iom*Z(iw)+gl(2,2)) ) /iom ;
      
    end
  
  end % ... of frequency loop
  
  rk = unwrap(real(k));
  ik = imag(k);
  
  [km,ikm] = max(rk);
  wmax_num(p) = w(ikm);
  
    %------------------
    % Dispersion curves
  
  if FIGS(1)
    figure(FIGS(1))
    subplot(ceil(np/2),2,p)
    plot(rk,w,ik,w,'--',[0 km],[0 km],'k:', km,w(ikm),'ro',km,w(ikm)/2,'r+')
    if p==np | p==np-1, xlabel('Re(k h) and Im(k h)'), end
    if any([1:2:np]==p), ylabel('\omega h/c'), end
    if p==1, legend('Re(k h)','Im(k h)','exact',4), end
    
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    text(xl(2)/2,0.87*yl(2),sprintf('NGLL = %u',NGLL),'HorizontalAlignment','center')
    
    % two scales in the same plot
    v = get(gca,'position');
    npwl = [10 6 5 4 3 2];
    box off
    axes('position',v,'Ytick',[],'Xlim',xl/(2*pi*P),...
         'xtick',1./npwl,'xticklabel',npwl,...
         'tickdir','in','XAxisLocation','top','color','none')
    if p==1 | p==2, xlabel('Nodes per wavelength'), end
    axes('position',v,'ytick',[],'xtick',[],'color','none','box','on')
  end
  
    %------------------
    % Impedance curves
  if FIGS(3)
    figure(FIGS(3))
    subplot(ceil(np/2),2,p)
    loglog(w,abs(Z))
    if p==np | p==np-1, xlabel('\omega h/c'), end
    if any([1:2:np]==p), ylabel('Z = T/V bottom'), end
    xl = get(gca,'xlim'); 
    xl = log10(xl);
    xl = 10^(xl(1)+0.5*(xl(2)-xl(1)));
    yl = get(gca,'ylim'); 
    yl = log10(yl);
    yl = 10^(0.6*yl(2));
    text(xl,yl,sprintf('NGLL = %u',NGLL),'HorizontalAlignment','center')
  end
  
end % ... of P loop

%------------------
% Maximum frequency stuff ==> stability
if FIGS(2)
  figure(FIGS(2))
  clf
  if KIND=='gll'
    subplot(121)
    plot(Ps+1,dxm_num,'o',Ps+1,dxm,'--')
    xlabel('NGLL')
    ylabel('\Delta x_{min} /h')
    legend('exact','4 /(NGLL^2-1)',1)
    subplot(122)
    plot(Ps+1,wmax_num,'o',Ps+1,wmax,':',Ps+1,7/3./dxm_num,'--')
    xlabel('NGLL')
    ylabel('\omega_{max} h/c')
    legend('exact','7/12 (NGLL^2-1)','7/3 h/\Delta x_{min}',2)
  
  else
    subplot(121)
    plot(Ps+1,dxm_num./dxm,'-o')
    xlabel('NGLL')
    ylabel('\Delta x_{min}^{prolate} /\Delta x_{min}^{SEM} ')
    subplot(122)
    plot(Ps+1,wmax_num./wmax,'-o')
    xlabel('NGLL')
    ylabel('\omega_{max}^{prolate} / \omega_{max}^{SEM}')
  
  end
end
