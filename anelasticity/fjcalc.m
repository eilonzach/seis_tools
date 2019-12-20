function [qinv,gg,ks]=fjcalc(T0, d, f, P)
%function [qinv,gg,ks]=fjcalc(T0, d, f, P)
% Calculate G, Q, ... vs. temperature from Jackson & Faul 2005
%   (EPSL 234, 119-134)
% Input:
%   T0 = temperature (C)
%   d = grain size (meters)
%   f = frequency (Hz)
%   P = pressure (GPa)    
% Output:  qinv   gg  kk
%  qinv = 1/Q
%  gg = shear modulus (GPa)
%  ks = adiabatic bulk modulus (GPa) 


w=2*pi*f;
T=T0+273;

% Their parameters in table 1
dr=1.e-5;    % m
Tr=1223;     % K
Ju=0.0149;    % 1/GPa
dlJudt=9.1e-4;  %1/K
Delta=1.4;
alph=.27;
tauLR=3.981e-3;
tauHR=5.26E6;
tauMR=4.31e6;
mA=1.09;
mJ=0.16;
mV=2.1;
E=5.05E5;    % J/mol
Vstar=1.2e-5;  %m^3/mol
dGdT0=-.0136;   %GPa/K
dGdP0=1.8;
R=8.314472;       % gas const., J/K/mol

% MY bulk modulus (anharmonic variations only)
dKsdT = -0.0204;   % GPa/K    for FO90 olivine
Ks0 = 123.;     %GPa, at 3 GPa pressure and 1200C
dKsdP = 5.03;    % dKs/DP, same conditions
ks = Ks0 + (T0-1200).*dKsdT + (P-3).*dKsdP;

% Pressure?
Jup=1./(1./Ju + P.*dGdP0);
if (T < Tr)
    dlnJu= -dGdT0.*(T - Tr).*Jup;   % low-T anharmonic only; ; see email from Uli 7/18/07)
else
  % viscoelastic and anelastic temperature effects
    dlnJu=dlJudt.*(T-Tr)./( (d/dr).^mJ);
end
expfac=exp( (E/R).*(1./T - 1./Tr) + (P.*1E9).*Vstar./(R.*T));
tauL=tauLR.*(d./dr).^mA.*expfac;
tauH=tauHR.*(d./dr).^mA.*expfac;
tauM=tauMR.*(d./dr).^mV.*expfac;

taufac=alph.*Delta./(tauH.^alph - tauL.^alph);
% try iterative quadrature   (still seems to be a problem with high freq
%    w*tauH > 1 say)
FINT1 = @(x) (x.^(alph-1))./(1+(w.*x).^2);
int1=taufac.*quadl(FINT1, tauL, tauH);

FINT2 = @(x) (x.^alph)./(1+(w.*x).^2);
int2=w.*taufac.*quadl(FINT2, tauL, tauH);
J1f=1 + dlnJu + int1;
J2f=int2 + 1./(w.*tauM);
qinv=J2f/J1f;
gg=1./Jup./sqrt(J1f.*J1f + J2f.*J2f);



return
