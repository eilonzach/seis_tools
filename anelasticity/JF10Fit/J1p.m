function ip1 = J1p(tau,omega,tauP)
global sig

ip1 = (1./tau).*exp(-0.5*(log(tau./tauP)./sig).^2)./(1+(omega.*tau).^2);