%% Simple example illustrating usage of JF 10 codes
T = 1200; % temperature in C
Z = 50;  % depth in km
frq = 1; % frequency in Hz
gs = 0.001; % grain size in m
% Probably don't change:
vfac = 1; % modification to viscosity prefactor - for melt/water 
% computed on the basis of the above, but you can sub in own values
P = Z/32; % pressure in GPA
omega = 2*pi*frq; 


%% Anharmonic velocities
% should really calculate these from P,T...
Vp_anh = 8.2e3; % Compressional velocity, m/s
Vs_anh = 4.3e3; % Shear velocity, m/s
rho = 3.3e3; % density, kg/m^3
G = Vs_anh.^2*rho; % elastic shear modulus, Pa
K = Vp_anh.^2*rho - (4/3)*G; % elastic bulk modulus, Pa

%% Using J&F 2010
[J1,J2]=creep10_GA(T+273,gs,P, omega*ones(size(P)),vfac); 




%% Results
qinv = J2./J1; % inverse Q
gg=G./sqrt(J1.^2 + J2.^2); % anelastic shear modulus

% Calculate
Qs = 1./qinv;
Qp = (9/4)*Qs; % using classic relationship
Vs = sqrt(gg./rho);
Vp = sqrt((K + 1.333*gg)./rho);

% output
fprintf('For T = %4.0f C, P = %4.2f GPa, gs = %4.2f mm, %4.2f Hz\n',...
    T,P/1e9,gs*1e3,frq)
fprintf('  Vs = %5.3f km/s   (reduced from anh value of %5.3f km/s)\n',...
    Vs/1e3,Vs_anh/1e3)
fprintf('  Vp = %5.3f km/s   (reduced from anh value of %5.3f km/s)\n',...
    Vp/1e3,Vp_anh/1e3)
fprintf('  Qs = %5.1f \n',...
    Qs)





