function int1 = J1anel(tau,omega)
global alpha
int1 = tau.^(alpha-1)./(ones(1,length(tau))+omega^2*tau.^2);
