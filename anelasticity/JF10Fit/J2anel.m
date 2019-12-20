function int2 = J2anel(tau,omega)
global alpha
int2 = tau.^(alpha)./(ones(1,length(tau))+omega^2*tau.^2);