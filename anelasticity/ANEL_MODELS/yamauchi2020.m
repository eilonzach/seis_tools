function [ J1,J2,tauM,qinv,G ] = yamauchi2020( T,d,P,omega,vfac,Gu,phi)
% [ J1,J2,tauM,qinv,G ] = takei2017( T,d,P,omega,vfac,Gu,phi)
% Inputs:
%   T = temperature (K)
%   d = grain size (meters)
%   P = pressure (GPa)    
%   omega = angular frequency (rad/s)
%   vfac = viscosity prefactor
%   Gu = unrelaxed shear modulus (GPa) if you want a non-default one. 
%   phi = volumetric melt fraction (porosity)
% Output: 
%  J1 = storage compliance (1/Pa)
%  J2 = loss  compliance (1/Pa)
%  tauM = Maxwell time (s)
%  qinv = inverse Q, 1/Q = J2/J1
%  G  = shear modulus (Pa)
% 
%%%%%
% pars = parameter struct

% this function derived from Yamauchi, Hatsuki, and Yasuko Takei.
% "Application of a Premelting Model to the Lithosphere‐Asthenosphere
% Boundary." Geochemistry, Geophysics, Geosystems 21, no. 11 (2020):
% e2020GC009338.
% 
% Z.E. 04/2022 modified 


if nargin<7, phi=0; end

d_mm = d*1000;

% Table 1
pars.AB         = 0.664;
pars.alphaB     = 0.38;
pars.tauP       = 6e-5;
pars.gamma      = 5;
pars.Teta       = 0.94;
pars.lambda     = 0;

% Table 2
pars.muU0       = 72.45;      % GPa unrelaxed compliance
pars.dmuU_dT    = -1.094e-2;   % GPa/K compliance variation with T
pars.dmuU_dP    = 1.987;      % dim'less compliance variation with P
pars.H          = 462.5e3;    % J/mol activation energy
pars.R          = 8.314;      % J/mol/K universal gas constant
pars.Omega      = 7.913e-6;   % m^3/mol activation volume
pars.dTm_dz     = 1.018;      % K/km solidus slope
pars.Tm_50km    = 1326 + 273; % K melting temperature
pars.km_per_GPa = 29.94;      % km/GPa pressure--depth conversion 
pars.log10_etar = 21.794;     % Pa-sec reference viscosity
pars.Tr         = 1200 + 273; % K reference temperature
pars.Pr         = 1.5;        % GPa reference pressure
pars.dr         = 0.5;        % mm reference grain size
pars.dm         = 1;          % dim'less grain-size exponent

% solidus & homologous temperature
depth_km = P*pars.km_per_GPa;
Tm_K = pars.Tm_50km + pars.dTm_dz*(depth_km - 50);
Th = T./Tm_K;

% compute viscosity reduction just below solidus, eqn. (29)
Aeta = 1;
if Th >= pars.Teta
  if Th < 1 % below melting point
      Aeta = exp(-(Th-pars.Teta)./(Th*(1-pars.Teta))*log(pars.gamma));
  else % DO include a melting effect, in line with Mei et al, 2002
      Aeta = exp(-pars.lambda*phi)/pars.gamma; 
  end
end
Aeta = Aeta.*vfac; % cheating adding an extra factor

% compute viscosity, eqn. (28)
eta = 10^pars.log10_etar.*(d_mm/pars.dr).^pars.dm.* ...
    exp(pars.H/pars.R.*(1./T - 1/pars.Tr)).* ...
    exp(pars.Omega/pars.R*(P./T - pars.Pr/pars.Tr)*1e9).* ...
    Aeta;

% compute shear modulus in GPa, eqn (27)
% reference conditions are 0 degC, 0 Pa
if nargin<6 || isempty(Gu)
  Gu = pars.muU0 + pars.dmuU_dT*(T-273) + pars.dmuU_dP*(P-0);
end

Ju = 1./(Gu*1e9); % Pa

% compute Maxwell time, seconds
tauM = Ju.*eta;

% compute peak width, eqn (25)
sigmaP = 4;
if Th>=0.92
  if Th>=1
      sigmaP = 7;
  else
      sigmaP = 4+37.5*(Th-0.92); 
  end
end

% compute peak amplitude, eqn (24)
AP = 0.01;
if Th>=0.91
  if Th>=1
      AP = 0.03; % + beta(phi)
  elseif Th>=0.96
      AP = 0.03;
  else
      AP = 0.01 + 0.4*(Th-0.91); 
  end
end

% compute normalised period, eqn (21)
p = 1./(omega.*tauM);

% compute the real and imaginary parts of the compliance, eqn (26)
J1 = Ju.*(1 + pars.AB*p.^pars.alphaB/pars.alphaB + ...
       sqrt(2*pi)/2*AP*sigmaP*(1 - erf( log(pars.tauP./p)/(sqrt(2)*sigmaP) )));
J2 = Ju.*(pi/2*(pars.AB*p.^pars.alphaB + ...
             AP*exp(-(log(pars.tauP./p)).^2/(2*sigmaP^2))) + p);

% form complex compliance
J = J1 + 1i*J2;

qinv=J2./J1;
G=abs(1./J)./1e9; % GPa
  