function [ K,G,M ] = moduli_from_VpVsRho( Vp,Vs,rho )
%[ K,G,M ] = moduli_from_VpVsRho( Vp,Vs,rho )

G = rho.*Vs.^2;
M = rho.*Vp.^2;
K = M - (4/3)*G;



end

