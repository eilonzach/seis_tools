function [ J1,J2,tauM,qinv,G ] = PM2013( T,d,P,omega,vfac,Gu )
% [ J1,J2,tauM,qinv,G ] = PM2013( T,d,P,omega,vfac,Gu )
%
%   Application of the anelastic scaling law of Priestley and McKenzie 2013
% 
% INPUTS:
%  T       temperature, in Kelvin
%  d       grain size, in m
%  P       pressure, in GPa
%  omega   angular frequency, in rad/s
%  vfac    viscosity reduction factor [optional - 1 if not given]
%  Gu      anharmonic shear modulus, in GPa [optional - calc. if not given]
% OUTPUTS:  
%  J1 = storage compliance (1/Pa)
%  J2 = loss  compliance (1/Pa)
%  tau_ss = steady-state relaxation time (s) = (viscos_ss)/Gu
%  qinv = inverse Q, 1/Q = J2/J1
%  G  = shear modulus (GPa)
% 
%  Z.E. 05/2017

%% parameters
mu0     = 72.66;        % GPa
dmudT   = -0.00871;     % GPa/K
dmudP   = 2.04e-9;         % GPa/Pa
Eta0    = 10^(22.38);   % Pa s
Ea      = 420.9e3;      % J/mol
Va      = 7.81e-6;      % m^3/mol

Pr      = 1.5e9;        % Pa
Tr      = 1473;         % K

dr      = 0.004;        % m (inferred reference grainsize at which Eta0 is implicitly calculated

R       = 8.3144621;    % ideal gas constant

%% Pre-things
if nargin<4 || isempty(omega)
    omega = 2*pi;
end
if any(d==0) || isempty(d)
    d = dr;    % m
end
if nargin<5 || isempty(vfac)
    vfac = 1;
end
if nargin<7 || isempty(rho)
    rho = 3300;     % kg/m^3
end

f = omega/2/pi;
P = P(:)*1e9;
T = T(:);

    
%% Anharmonic modulus/compliance
if nargin<6 || isempty(Gu)
Gu = mu0 + dmudT*T + dmudP*P;
end
Gu = Gu*1e9; % in Pa
Ju = 1./Gu;

%% viscosity
astar = exp((Ea + Pr*Va)/(R*Tr))./exp((Ea + P*Va)./(R*T));
Eta = Eta0./astar;

%% account for grain size, vfac
Eta = Eta.*vfac.*(d/dr).^3;

%% Maxwell time (diffusion creep)
tauM = Ju.*Eta;
fn = tauM.*f;
taun = 1./(2*pi*fn);

%% Complex compliances
J1 = Ju.*(1./F_master_M11(fn));
J2 = Ju.*((pi/2)*Xn_master_M11(taun) + taun );

J = (J1 + 1i*J2);

%% anelastic shear mod, atten
qinv = J2./J1;
G=1/sqrt(J1.*J1 + J2.*J2);



end




function F = F_master_M11(fn)
    a0 =  0.55097;
    a1 =  0.054332;
    a2 = -0.0023615;
    a3 = -5.7175e-5;
    a4 =  9.9473e-6;
    a5 = -3.4761e-7;
    a6 =  3.9461e-9;
    fr =  1e13;
    lnf = log(fn);
    
    F = a0 + a1.*lnf + a2.*lnf.^2 + a3.*lnf.^3 + a4.*lnf.^4 + a5.*lnf.^5 + a6.*lnf.^6;
    F(fn>fr) = 1;
end

function Xn = Xn_master_M11(tau)
    Xn = nan(size(tau));
    hitau = tau>=1e-11;
    lotau = tau<1e-11;
    
    y = 0.39 - 0.28./(1 + 2.6*(tau.^0.1));
    
    Xn(hitau) = 0.32.*tau(hitau).^(y(hitau));
    Xn(lotau) = 1853.*sqrt(tau(lotau));
end

    
    