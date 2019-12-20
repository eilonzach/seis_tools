function ip2 = J2p(tau,omega,tauP)
global sig

ip2 = exp(-0.5*(log(tau./tauP)./sig).^2)./(1+(omega.*tau).^2);