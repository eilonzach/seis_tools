function [J1,J2] = creep10mod(Te,gs,pres,omega,dvis)
%function [J1,J2] = creep10mod(Te,gs,pres,omega,dvis)
%calculate real and imaginary parts of creepfunction, eq. 7 JF10
%params from Jackson and Faul, PEPI, 2010, Table 2, note p
% from U Faul 
%   Ts =   Temperature  (K)
%   gs =   grain size  (m)
%   pres = Pressure  (GPa)
%   omega = angular frequency (1/s)
%   [dvis] = added correction factor to viscosity/relaxation time, e.g. (COH/Cr)^-r etc
%           (GA) assumed same for all elements, as is E*, V*
%  Ts, omega should be matching vectors (or all are)
% OUTPUT:  J1, J2 are real, imaginary compliance, normalized to unrelaxed Ju

global alpha sig% alpha is used in J1anel and J2anel, sig in J1p and J2p
Tr = 1173; iTr=1/Tr; %reference temperature in K
Pr = 0.2; PT = Pr/Tr; %reference pressure in GPa
gsr = 1.34E-5; % reference grain size in m
deltaB = 1.04; % background relaxation strength,
alpha = 0.274; % background frequency exponent
%reference values for relaxation times
tauLo = 1E-3; tauHo = 1E7; tauMo = 3.02E7; 
ma = 1.31;  % anelastic grain size exponent
mv = 3;     % viscous grain size exponent
EB = 3.6E5; % activation energy for background and peak
AV = 1E-5;  % activation volume
R = 8.314; 
AVR = AV/R; ER = EB/R; gr = gs/gsr;

% peak parameters:
tauPo = 3.98E-4; % reference peak relaxation time,
deltaP = 0.057;  % peak relaxation strength pref
sig = 4;         % peak width
cp = deltaP * (2*pi)^(-0.5)/sig; %peak integration const.

% relaxation times eqs. 9 and 10
% GAA:  Added the ".*1.e9" to fix a bug in original code  11/13
taut = dvis.*exp((ER)*(1./Te-iTr)).*exp(AVR*((pres./Te)-PT).*1.e9);
tauH = tauHo*gr.^ma.*taut;
tauL = tauLo*gr.^ma.*taut;
tauP = tauPo*gr.^ma.*taut;
tauM = tauMo*gr.^mv.*taut;
%initialize arrays
ij1 = zeros(1,length(Te)); ij2 = zeros(1,length(Te)); on = ones(size(Te));
ip1 = zeros(1,length(Te)); ip2 = zeros(1,length(Te));

%integration for peak wrt to dtau 0 - inf;
for ii = 1 : length(Te)
  ij1(ii) = quadl(@(tau)J1anel(tau,omega(ii)),tauL(ii),tauH(ii));
  ij2(ii) = quadl(@(tau)J2anel(tau,omega(ii)),tauL(ii),tauH(ii));
  ip1(ii) = quadgk(@(tau)J1p(tau,omega(ii),tauP(ii)),0,inf);
  ip2(ii) = quadgk(@(tau)J2p(tau,omega(ii),tauP(ii)),0,inf);
end

Jb1 = alpha*deltaB*ij1./(tauH.^alpha-tauL.^alpha);
Jb2 = omega*alpha*deltaB.*ij2./(tauH.^alpha-tauL.^alpha);
Jp1 = cp.*ip1; 
Jp2 = cp.*omega.*ip2;

J1 = on + Jb1 + Jp1;
J2 = (Jb2 + Jp2) + 1./(omega.*tauM);

return