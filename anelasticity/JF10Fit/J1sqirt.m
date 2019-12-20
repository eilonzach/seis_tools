function isq1 = J1sqrt(tau,omega,tauSQ)
global sigsq

isq1 = (1./tau).*exp(-0.5*(log(tau./tauSQ)./sigsq).^2)./(1+(omega.*tau).^2);