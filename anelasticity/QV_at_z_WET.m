function [ Qs,Qp,Vs,Vp,gs ] = QV_at_z_WET( age,Z,COH,frq,rate)
% [ Qs,Qp ] = Qs_z( age )
%   Function to calculate a Q(z) profile, assuming a certain geotherm and
%   then using the Jackson and Faul, 2010 treatment to predict Qs etc.

if nargin < 2
    Z = [10:5:200]';
end
if nargin < 3
    COH = 50;       % PPM COH in H/Si
end
if nargin < 4
    frq = 1;
end
if nargin < 5
    rate = 5e-14;      % strain rate, s^-1
end

omega = 2*pi*frq;

%% parms
% Parameters for grain size/stress calculation:  strain rate, P, COH, T  array
% dealing with water
COH = COH;       % up to 75*16;  PPM COH in H/Si
COHref = 50;     % reference background ppm H/Si
rh     = 1.0;     % power law from Hirth & Kohlstedt 2003 for Coh
dvis_OH=(COH/COHref)^(-rh);

addpath('/Users/zeilon/Documents/MATLAB/geoff/VpVsQ')
addpath('/Users/zeilon/Documents/MATLAB/geoff/VpVsQ/FinalScripts-AbersQ2014-2/JF10Fit')

TempZ = geotherm( age,'plate',Z,1350);
P = Z/32;

%% Get anharmonic moduli from HeFESTo
% % Geoff's bulk modulus (anharmonic variations only)
% dKsdT = -0.0204;   % GPa/K    for FO90 olivine
% Ks0 = 123.;     %GPa, at 3 GPa pressure and 1200C
% dKsdP = 5.03;    % dKs/DP, same conditions
% ks = Ks0 + (TempZ-1200).*dKsdT + (P-3).*dKsdP;

[rho,ks,G,VSh,VSv,VPh,VPv] = HeFESTo_getval( TempZ, P );


%% calc anelasticity
[~,gs]=grainsizeae(TempZ,rate,P*1e9,COH);     %TS=stress, d=grain size in m vs. T
[J1,J2]=creep10_GA(TempZ+273,gs,P, omega*ones(size(P)),dvis_OH); 
Qs = J1./J2;

gg=G./sqrt(J1.^2 + J2.^2);



Qp = (9/4)*Qs; % using classic relationship
Vs = sqrt(gg*1e9./rho);
Vp = sqrt((ks + 1.333*gg).*1e9./rho);

% figure(99), hold on
% plot(Qs,-Z)
% plot(Vs/50,-Z)
% plot(Vp/50,-Z)
% xlim([0 500])



end

