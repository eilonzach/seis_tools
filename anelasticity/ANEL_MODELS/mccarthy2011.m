function [J1,J2,tauM,qinv,G]=mccarthy2011(T, d, P, omega, vfac, Gu)
% (1) function [qinv,gg]=atten_mccarth11(T0, d, f, P)
% (2) function [qinv, gg]=atten_mccarth11(visc,Gun,f);
% Calculate G, Q, ... vs. temperature etc from McCarthy et al. 2011
%   (JGR)
%  GAA 12/11
% INPUT (1):  (Used for the Olivine flow law extrap. in their paper)
%   T = temperature (K)
%   d = grain size (meters)
%   omega = angular frequency (rad/s)
%   P = pressure (GPa)    
%   vfac = viscosity prefactor
%   Gu = unrelaxed shear modulus (GPa)
% INPUT (2):  (user-provided viscosity/elasticity, same master curve)
%   visc = viscosity array, Pa s
%   Gun = unrelaxed shear modulus, GPa
%   omega = angular frequency (rad/s)
%      to get Gun for olivine (Isaak), try Gun=82 - 0.0136*(T+273) + 1.8*P  
% OUTPUTS:  
%  J1 = storage compliance (1/Pa)
%  J2 = loss  compliance (1/Pa)
%  tau_ss = steady-state relaxation time (s) = (viscos_ss)/Gu
%  qinv = inverse Q, 1/Q = J2/J1
%  G  = shear modulus (GPa)
% 
%  Z.E. 05/2017

f = omega/2/pi;

% Their parameters:  reference state for olivine after Faul/Jackson
dr=1.e-3;    % grain size, m
m=3;         % grain size exponent ??!! (this is theory; not clear is consistent)
Tr=1473;     % K
Pr=0.0;       % GPa
E=5.05E5;    % J/mol
V=1.20e-5;  %m^3/mol
visc0=6.6E19;   % ref viscosity, Pa s  
R=8.314472;       % gas const., J/K/mol
%anharmonic params for olivine from Isack/Anderson
Gur = 82;    % GPa  -- seems stiff vs. 66 in FJ05/JF10
dGdT0=-.0136;   %GPa/K
dGdP0=1.8;
%Anelastic: polynomials for Ju/J1  (eqn 26)
a=[.55097 .054332 -.0023615 -5.7175e-5  9.9476e-6 -3.4761e-7 3.9461e-9];
tnmin=1.e-11;    % min period below which elastic
%tnmin=1.e-13;    % fake to make go lower

% Unrelaxed Pressure/etc deriv's  -- GPa
if nargin<6 || isempty(Gu)
    Gu = Gur+ P.*dGdP0 + T.*dGdT0;
end
Ju=1./Gu;
% Reference maxwell time
tauM=vfac*(Ju.*1.E-9).*visc0.*((d./dr).^m).*...
      exp(E./R.*(1./T - 1./Tr)).*exp(V./R.*(P./T-Pr./Tr).*1.e9);
% alternate input format
if (nargin()<4)
    Ju=1./d;
    visc=T;
    tauM=(Ju.*1.E-9).*visc;
end
% Master variable
fn=f.*tauM;   % normalized frequency;
tn=1./(2.*pi.*fn);     % normalized period
Xn = 0.32.*tn.^(0.39-0.28./(1+2.6.*tn.^0.1));    %eqn. 25
k11=find(tn<tnmin);
if(~isempty(k11)), Xn(k11)=1853.*tn(k11).^0.5; end;
% Ju/J1  from eqn. 26
juj1=zeros(size(fn));
for k=0:6
  juj1=juj1+a(k+1).* (log(fn).^k);
end
k13=find(fn>1.E13);
if (~isempty(k13)), juj1(k13)=1; end,

J1f=1./juj1;
J2f=0.5.*pi.*Xn + tn;

J1 = Ju*J1f;
J2 = Ju*J2f;

% rest is standard

qinv=J2./J1;
G=1./sqrt(J1.*J1 + J2.*J2);

return
