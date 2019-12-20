function [J1,J2,tau_ss] = creep10_GA(Te,gs,pres,omega,vfac, ipk)
%function [J1,J2] = creep10_GA(Te,gs,pres,omega,vfac, ipk)
%function [J1,J2,tau_ss] = creep10_GA(Te,gs,pres,omega,vfac, ipk)
%calculate real and imaginary parts of creepfunction, eq. 7 JF10
%params from Jackson and Faul, PEPI, 2010, Table 2, note p
% Te = temperature (K)
% gs = grain size (m)
% pres = pressure (GPa)
% omega = angular frequency (rad/s)
%  This version allows vfac = scale factor to all relaxation times
% Te, gs, pres, omega,vfac should all be matching arrays/same length
% ipk=1 to turn off the HT Peak THIS IS ALWAYS ON/KEPT FOR COMPATIBILITY
% J1, J2 are real, imaginary parts of complex compliance, normalized to Gu
%  tau_ss = steady-state relaxation time, s = (viscos_ss)/Gu

IDOPEAK=0;
if (nargin()>5)
    if ipk>0, IDOPEAK=1; end;
end
if numel(omega)==1
    omega = omega*ones(size(Te));
end
if numel(gs)==1
    gs = gs*ones(size(Te));
end

global alpha sig; % alpha is used in J1anel and J2anel, sig in J1p and J2p
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

% peak parameters:
tauPo = 3.98E-4; % reference peak relaxation time,
deltaP = 0.057;  % peak relaxation strength pref
sig = 4;         % peak width

% TEST TEST TEST - this is Multiple Sol Gel only (line m Table 2 JF10)
%  this line:  fixes mV=3 (creep d^mV) and fits peak
%ma = 0.887; deltaB=1.91; alpha=.301; tauMo=1.32E6; EB=3.48E5; deltaP=.072; tauPo=1.995E-4;


AVR = AV/R; ER = EB/R; gr = gs/gsr;
cp = deltaP * (2*pi)^(-0.5)/sig; %peak integration const.

% relaxation times eqs. 9 and 10
% GAA:  Added the ".*1.e9" to fix a bug in original code  11/13
% GAA - add vfac to original, uniform global scaling
taut = vfac.*exp((ER)*(1./Te-iTr)).*exp(AVR*((pres./Te)-PT).*1.e9);
tauH = tauHo*gr.^ma.*taut;
tauL = tauLo*gr.^ma.*taut;
tauP = tauPo*gr.^ma.*taut;
tauM = tauMo*gr.^mv.*taut;
%initialize arrays
ij1 = zeros(size(Te)); 
ij2 = zeros(size(Te)); 
on = ones(size(Te));
ip1 = zeros(size(Te)); 
ip2 = zeros(size(Te)); 
%integration for peak wrt to dtau 0 - inf;
for ii = 1 : length(Te)
  ij1(ii) = quadl(@(tau)J1anel(tau,omega(ii)),tauL(ii),tauH(ii));
  ij2(ii) = quadl(@(tau)J2anel(tau,omega(ii)),tauL(ii),tauH(ii));
  ip1(ii) = quadgk(@(tau)J1p(tau,omega(ii),tauP(ii)),0,inf);
  ip2(ii) = quadgk(@(tau)J2p(tau,omega(ii),tauP(ii)),0,inf);
end

Jb1 = alpha*deltaB*ij1./(tauH.^alpha-tauL.^alpha);
Jb2 = omega.*alpha*deltaB.*ij2./(tauH.^alpha-tauL.^alpha);
Jp1 = cp.*ip1; 
Jp2 = cp.*omega.*ip2;

J1 = on + Jb1 + Jp1;
J2 = (Jb2 + Jp2) + 1./(omega.*tauM);

if (nargout()>2)
    tau_ss=tauM;
end

return