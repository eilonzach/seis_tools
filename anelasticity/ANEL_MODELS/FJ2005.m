function [J1,J2,tauM,qinv,G]=FJ2005(T, d, P, omega, vfac,Gu)
%[J1,J2,tauM]=jf2005(T, d, P, omega, vfac)
% Calculate G, Q, ... vs. temperature from Jackson & Faul 2005
%   (EPSL 234, 119-134)
% Inputs:
%   T = temperature (K)
%   d = grain size (meters)
%   P = pressure (GPa)    
%   omega = angular frequency (rad/s)
%   vfac = viscosity prefactor
% Output:  qinv   gg  kk
%  J1 = storage compliance (1/Pa)
%  J2 = loss  compliance (1/Pa)
%  tauM = Maxwell time (s)
%  qinv = inverse Q, 1/Q = J2/J1
%  G  = shear modulus (GPa)
% 
%  Z.E. 05/2017
%  modified to include an alteration to the pseudoperiod prefactor to
%  account for viscosity change

f = omega./2./pi;

w=2*pi*f;
T0 = T-273;

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

% Unrelaxed compliance
if nargin<6 || isempty(Gu)
    Ju=1./(1./Ju + P.*dGdP0);
else
    Ju = 1./Gu;
end

expfac=vfac.*exp( (E/R).*(1./T - 1./Tr) + (P.*1E9).*Vstar./(R.*T));
tauL=tauLR.*(d./dr).^mA.*expfac;
tauH=tauHR.*(d./dr).^mA.*expfac;
tauM=tauMR.*(d./dr).^mV.*expfac;

taufac=alph.*Delta./(tauH.^alph - tauL.^alph);
% try iterative quadrature   (still seems to be a problem with high freq
%    w*tauH > 1 say)
FINT1 = @(x) (x.^(alph-1))./(1+(w.*x).^2);
int1 = zeros(length(taufac),1);
for ii = 1:length(taufac)
    int1(ii)=taufac(ii).*quadl(FINT1, tauL(ii), tauH(ii));
end
    

FINT2 = @(x) (x.^alph)./(1+(w.*x).^2);
int2 = zeros(length(taufac),1);
for ii = 1:length(taufac)
    int2(ii)=w.*taufac(ii).*quadl(FINT2, tauL(ii), tauH(ii));
end

% already accounted for anharmonic effects
% if (T < Tr)
%     dlnJu= -dGdT0.*(T - Tr).*Jup;   % low-T anharmonic only; ; see email from Uli 7/18/07)
% else
%   % viscoelastic and anelastic temperature effects
%     dlnJu=dlJudt.*(T-Tr)./( (d/dr).^mJ);
% end
% J1 = 1 + dlnJu + int1;  J1 = abs(J1);

J1f = 1 + int1;  J1f = abs(J1f);
J2f = int2 + 1./(w.*tauM);

J1 = Ju*J1f;
J2 = Ju*J2f;

% rest is standard

qinv=J2./J1;
G=1./sqrt(J1.*J1 + J2.*J2);


return
