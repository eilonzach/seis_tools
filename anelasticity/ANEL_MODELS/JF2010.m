function [J1,J2,tauM,qinv,G] = JF2010(T,d,P,omega,vfac, Gu)
%[J1,J2,tauM_ss,qinv,G] = JF2010(T,d,P,omega,vfac, Gu)
%calculate real and imaginary parts of creepfunction, eq. 7 JF10
%params from Jackson and Faul, PEPI, 2010, Table 2, note p
% INPUTS:
%   T = temperature (K)
%   d = grain size (meters)
%   P = pressure (GPa)    
%   omega = angular frequency (rad/s)
%   vfac = viscosity prefactor
%   Gu = unrelaxed shear modulus (GPa)
% OUTPUTS:
%  J1 = storage compliance (1/Pa)
%  J2 = loss  compliance (1/Pa)
%  tauM = Maxwell relaxation time (s) = (viscos_ss)/Gu
%  qinv = inverse Q, 1/Q = J2/J1
%  G  = shear modulus (GPa)
% 
%  Z.E. 05/2017
%  This version allows vfac = scale factor to all relaxation times
% Te, gs, pres, omega,vfac should all be matching arrays/same length
% ipk=1 to turn off the HT Peak THIS IS ALWAYS ON/KEPT FOR COMPATIBILITY
% J1, J2 are real, imaginary parts of complex compliance, normalized to Gu
%  tau_ss = steady-state relaxation time, s = (viscos_ss)/Gu

addpath('/Users/zeilon/Dropbox/MATLAB/seis_tools/anelasticity/JF10Fit');


if numel(omega)==1
    omega = omega*ones(size(T));
end
if numel(d)==1
    d = d*ones(size(T));
end
if nargin<6
    Gu=1;
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


AVR = AV/R; ER = EB/R; gr = d/gsr;
cp = deltaP * (2*pi)^(-0.5)/sig; %peak integration const.

% relaxation times eqs. 9 and 10
% GAA:  Added the ".*1.e9" to fix a bug in original code  11/13
% GAA - add vfac to original, uniform global scaling
taut = vfac.*exp((ER)*(1./T-iTr)).*exp(AVR*((P./T)-PT).*1.e9);
tauH = tauHo*gr.^ma.*taut;
tauL = tauLo*gr.^ma.*taut;
tauP = tauPo*gr.^ma.*taut;
tauM = tauMo*gr.^mv.*taut;
%initialize arrays
ij1 = zeros(size(T)); 
ij2 = zeros(size(T)); 
on = ones(size(T));
ip1 = zeros(size(T)); 
ip2 = zeros(size(T)); 
%integration for peak wrt to dtau 0 - inf;
for ii = 1 : length(T)
  ij1(ii) = quadl(@(tau)J1anel(tau,omega(ii)),tauL(ii),tauH(ii));
  ij2(ii) = quadl(@(tau)J2anel(tau,omega(ii)),tauL(ii),tauH(ii));
  ip1(ii) = quadgk(@(tau)J1p(tau,omega(ii),tauP(ii)),0,inf);
  ip2(ii) = quadgk(@(tau)J2p(tau,omega(ii),tauP(ii)),0,inf);
end

Jb1 = alpha*deltaB*ij1./(tauH.^alpha-tauL.^alpha);
Jb2 = omega.*alpha*deltaB.*ij2./(tauH.^alpha-tauL.^alpha);
Jp1 = cp.*ip1; 
Jp2 = cp.*omega.*ip2;

J1f = on + Jb1 + Jp1;
J2f = (Jb2 + Jp2) + 1./(omega.*tauM);

J1 = J1f./Gu;
J2 = J2f./Gu;

% rest is standard

qinv=J2./J1;
G=1./sqrt(J1.*J1 + J2.*J2);

return