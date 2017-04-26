function [ pol,uP,uSV,uSH ] = pol_cmtpred(strike,dip,rake,M0,takeoff,foraz)
% [ pol,uP,uSV,uSH ] = pol_cmtpred(strike,dip,rake,takeoff,foraz)
%POL_CMTPRED Summary of this function goes here
% Predict the polarisation of an arriving P or S wave from the source
% mechanism parameters
%
% INPUTS (all in degrees)
% strike    -   strike of the fault plane
% dip       -   dip of the fault plane
% rake      -   rake of the slip on the fault plane
% takeoff   -   takeoff angle of the wave
% rayparam  -   ray parameter - depends on gcarc
% foraz     -   event to station azimuth
% ----------------------
% mu
% rho
% A
% a - alpha - P-wave velocity
% B - beta  - S-wave velocity
% r - distance from source to receiver (great circle? linear?)
% Tp        -   P travel time 
% Ts        -   S travel time

% %testing
% clc
% strike=10;
% dip=67;
% rake=-127;
% foraz=-90;
% takeoff=30;
%% Start here

strike = d2r(strike);
dip    = d2r(dip);
rake   = d2r(rake);
takeoff= d2r(takeoff);
foraz  = d2r(foraz);
u = 1;
R = 6371;
if ~exist('pspeed','var')==1, pspeed=5.8; end  % iasp91 surface p velocity is 5.8km/s
if ~exist('sspeed','var')==1, sspeed=3.36; end % iasp91 surface s velocity is 3.36km/s

fnormal2=[-sin(dip)*sin(strike);...
           sin(dip)*cos(strike);...
          -cos(dip)];

% wavespeeds
a = pspeed;     % P-wavespeed at surface 
B = sspeed;     % S-wavespeed at surface
% horizontal slowness
pa=sin(takeoff)/a;
pB=sin(takeoff)/B;
% vertical slowness
qa=sqrt(a^-2 - pa^2); % p-wave
qB=sqrt(B^-2 - pa^2); % s-wave


%% Ray directions
% slip, u
U = u*[cos(rake)*cos(strike) + cos(dip)*sin(rake)*sin(strike);...
       cos(rake)*sin(strike) - cos(dip)*sin(rake)*cos(strike);...
      -sin(rake)*sin(dip)];
% fault normal ( = nu)
V = [-sin(dip)*sin(strike);...
      sin(dip)*cos(strike);...
     -cos(dip)];
% P-wave direction ( = l = gamma)
dP = [sin(takeoff)*cos(foraz);...
      sin(takeoff)*sin(foraz);...
      cos(takeoff)];
% SV-wave direction ( = ?p)
dSV = [cos(takeoff)*cos(foraz);...
       cos(takeoff)*sin(foraz);...
      -sin(takeoff)];
% SH-wave direction ( = ?phi)
dSH = [-sin(foraz);...
       cos(foraz);...
       0];
%% Radiation patterns
Fp = 2*(dP'*V)*(dP'*U)/norm(U);

Fsv = ((dP'*V)*(U'*dSV) + (dP'*U)*(V'*dSV))/norm(U);

Fsh = ((dP'*V)*(U'*dSH) + (dP'*U)*(V'*dSH))/norm(U);

%% Far-field displacements
% factor on the outside - not relevant for the polarisation
% X = (mu*A)/(4*pi*rho*r)
dtP=1;
dtS=1;
uP  = Fp*norm(U)*(1/a^3)*dtP*dP;
uSV = Fsv*norm(U)*(1/B^3)*dtS*dSV;
uSH = Fsh*norm(U)*(1/a^3)*dtS*dSV; % N.B. is this a typo on the bottom?

% %% Radiation pattern for spherically symmetrical medium
% % factor on the outside - not relevant for the polarisation
% % Y = (mu*A)/(4*pi*sqrt(rho[source]*rho[receiver]*B(source)*B(receiver))*B[receiver[^2*R)
% dt=1; % properly should be t-Tp or Ts, but we don't care about this
uP = Fp*norm(U)*dP;         uP(abs(uP)<1e-10)=0;
uSV = Fsv*norm(U)*dSV;      uSV(abs(uSV)<1e-10)=0;
uSH = Fsh*norm(U)*dSH;      uSH(abs(uSH)<1e-10)=0;

%% Free surface correction
% Free surface coefficients
C1 = 2*B^-2 * (B^-2 - pa^2)/((B^-2 - pa^2)^2 + 4*pa^2 * qa * qB);
C2 = 4*B^-2 * qa * qB     /((B^-2 - pa^2)^2 + 4*pa^2 * qa * qB);
% surface displacement matrix
W = [-a*qa*C1  B*pa*C2  0
      a*pa*C2  B*qB*C1  0
         0        0     2];


%% Compute polarisation
uSV_surf = uSV(1:2)*W(2,2); % projection of uSV onto horizontal surface
uSH_surf = uSH(1:2)*W(3,3); % uSH at the horizontal surface
u_surf=uSV_surf+uSH_surf;
pol=r2d(atan2(u_surf(2),u_surf(1))); % this is the angle (in radians) from the SV direction


% % testing
% uP
% uSV
% uSH
% pol

end
