% Calculate Vs and Q vs depth for oceanic geotherms with Burgers model
% For parameters of experimental fit, see creep10.m
% 
clear all
% Parameters set for each run:
% distance from the ridge axis is equivalent to setting the age of occeanic
% lithosphere for the velocity of 5 cm below
dist = 5500000;
% dist = 500000 for age = 10Myr, 600000 12My, 1750000 35My, 5500000 110My

Tp = 1623; % mantle potential temperature for geotherm, K

%depth sampling, m
delz = [(5000:2000:197000),(200000:5000:400000)]; 
dk = -delz(2:length(delz))/1000; lz = length(delz);
% Period
% as a function of depth for surface waves(Forsyth, Geophys. Mon.,71,1992)
period = 3*delz/4200;
% fixed period for body waves
%period = ones(1,length(delz))*10;
omega = 2*pi./period;

fgs = 1E-3; % constant grain size in m
% for change of grain size with depth replace fgs with different gs below
% specify depth ranges and attempt to create smooth transitions
ugst = find(delz==165000); tgst = find(delz==350000);
gst(1:ugst)=ones(1,ugst)*fgs;
% make transition asymmetric with smaller initial slope
gst(ugst+1) = fgs;
gst(ugst+2:ugst+2+(lz-tgst))= ones(1,lz-(tgst-1))*fgs;
zg = [delz(1:ugst),200000, delz((tgst):lz)];
% smooth transition
gs = interp1(zg,gst,delz,'spline');
%plot(delz/1000,gs,'r')

GUR = 66.5; % anharmonic modulus at Tr and Pr, part of fit. 
% (Could be changed to simulate compositional variations.)
dGdT = -0.0136; dGdP = 1.8; %anharmonic pressure and temperature derivatives
% Bulk modulus
K0 = 129; dKdT = -0.016; dKdP = 4.6; % Bass 1995

Tr = 1173; Pr = 0.2; %reference temperature and pressure, part of fit
% parameters for T and density calculation (Turcotte and Schubert)
grav = 9.98; kappa = 1E-6; 
rho0 = 3310; 
betaa = 6E-12;

% adiabatic T and rho increase, Turcotte and Schubert, Ch. 4-16 and 4-27
Tad = Tead(delz,Tp);
% add conductive cooling
vel = 0.05/(60*60*24*365); age = (dist/vel)/(60*60*24*365*1E6);
Te = 273 + (Tad - 273) .* erf(delz/(2*sqrt(kappa*dist/vel)));
ff = find(Te+5 > Tad);

% calculate density and pressure
rr = 1 - rho0*grav*betaa*delz;
pres = -log(rr) / betaa;
presG = pres/1E9;
rho = rho0./rr;
GU = (GUR + (Te - Tr)* dGdT + (presG - Pr)* dGdP)*1E9;
Vs0 = sqrt(GUR*1E9/rho0);

int=find(Te>900); % calculate G and Q only where there are anelastic contributions
Ti=Te(int); gsi=gs(int); presi=pres(int); omegai=omega(int); di = -delz(int)/1000;
%Evaluate integrals as a function of temperature, grain size, pressure and
%frequency. J1 and J2 are the real and imaginary parts of the creep
%function.
[J1,J2] = creep10(Ti,gsi,presi,omegai);
Q = J1./J2;
Ga = (J1.^2+J2.^2).^(-0.5).*GU(int);
Gc = GU;
Gc(int) = Ga;

Vs = sqrt(Gc./rho);
Vsanh = sqrt(GU./rho);
   
KT = K0 + (dKdT*(Te-273));
KTp = KT + dKdP*presG;
Vp = sqrt((KTp*1E9+4*Gc/3)./rho);

figure(1)
plot(Vs(2:lz)/1000,dk,'r')%,Vsanh(2:lz)/1000,dk,'k')
xlabel('Vs, km/s')
ylabel('Depth, km')
axis([4.1 4.9 -400 0])
set(gcf,'paperposition',[1 4 9.3 7.5])
set(gcf,'papertype','usletter','paperunits','centimeter')

figure(2)
plot(Q,di,'r')
ylabel('Depth, km')
xlabel('Q')
axis([0 400 -400 0])
set(gcf,'paperposition',[1 4 5 4])
set(gcf,'papertype','usletter','paperunits','centimeter')
