function [ Qs,Qp,Vs,Vp ] = QV_at_z( age,Z,gs,frq, model,vfac )
% [ Qs,Qp ] = Qs_z( age )
%   Function to calculate a Q(z) profile, assuming a certain geotherm and
%   then using various models treatment to predict Qs etc.

if nargin < 2 || isempty(Z)
    Z = [5:5:200]';
else
    Z = Z(:);
end
if nargin < 3 || isempty(gs)
    gs = 0.001;
end
if nargin < 4 || isempty(frq)
    frq = 1;
end
if nargin < 5 || isempty(model)
    model = 'JF10';
end
if nargin < 6 || isempty(vfac)
    vfac = 1;
end

if numel(gs)==1
    gs = gs*ones(size(Z));
end

TempZ = geotherm( age,'halfspace',Z,1350);
P = Z/32;
omega = 2*pi*frq;

addpath('/Users/zeilon/Documents/MATLAB/geoff/VpVsQ')

%% Anharmonic from Hefesto
[rho,ks,G,VSh,VSv,VPh,VPv] = HeFESTo_getval( TempZ, P );

switch model
    case 'JF05'
        %% Using J&F 2005
        [J1,J2]=jf2005_ZE(TempZ+273,gs,P, omega*ones(size(P)),vfac); 
    case 'JF10'
        %% Using J&F 2010
        [J1,J2]=creep10_GA(TempZ+273,gs,P, omega*ones(size(P)),vfac); 
    case 'PM13'
        %% Using PM 2013 (adapted)
        [J1,J2]=PM13(TempZ+273,gs,P, omega*ones(size(P)),vfac); 
end



%% Results
qinv = J2./J1;
gg=G./sqrt(J1.^2 + J2.^2);

Qs = 1./qinv;
Qp = (9/4)*Qs; % using classic relationship
Vs = sqrt(gg*1e9./rho);
Vp = sqrt((ks + 1.333*gg)*1e9./rho);

% plot(Qs,-Z)
% xlim([0 500])



end

