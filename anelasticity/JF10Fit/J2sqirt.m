function isq2 = J2sqirt(tau,omega,tauSQ)
global sigsq

isq2 = exp(-0.5*(log(tau./tauSQ)./sigsq).^2)./(1+(omega.*tau).^2);