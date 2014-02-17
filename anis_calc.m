% function [ anis,phi ] = anis_calc( foraz,inc,pol,az,ap,anis )
%THIS SCRIPT MAY BE DEFUNCT
%   function to calculate effective anisotropy for shear waves propagating
%   through an anisotropic medium with hexagonal fabric described by an
%   azimuth and plunge of a symmetry axis, as well as a strength of
%   anisotropy, defined as (v1-v2)/(v1+v2) where v1 is the velocity of
%   particle motion in the symmetry axis vector and v2 is the velocity in
%   the two perpendicular directions. anis<0 corresponds to a slow symmetry
%   axis
%   ALL ANGLES IN DEGREES
% NOT READY FOR USE
% 
% INPUTS
%   foraz  = propagation azimuth of wave
%   inc    = inclination angle of wave (degrees from vertical)
%   pol    = polarisation of wave: angle from pure SV, between 0 and 180
%   phi    = azimuth of hexagonal symmetry axis 
%   ap     = plunge of hexagonal symmetry axis (degrees from horizontal)
%   anis   = percent (v1-v2)/(v1+v2) where v1 is the symmetry velocity
% 
% OUTPUTS
%  anis    = effective anisotropy, defined as (u1-u2)/(u1+u2) where u1 and
%               u2 are the maximum and minimum velocities in the plane
%               perpendicular to shear wave propagation direction
%  phi     = angle between the initial polarisation and the fast
%               polarisation
% 
% define some directions:
% 1) North
% 2) East
% 3) Vertical (down)
% 
% If az = 0, ap = 0, then the hexagonal symmetry axis is due N 

%% FOR TESTING
foraz  = 45;
inc    = 90;
pol    = 45; % = atand(SH/SV)
az     = 0;
ap     = 0;
anis   = 5;
%% END TESTING

anis = anis/100;


V = [1 0 0; 0 (1-anis)/(1+anis) 0; 0 0 (1-anis)/(1+anis)]; % velocities matrix
Rv1 = [cosd(az) sind(az) 0; -sind(az) cosd(az) 0; 0 0 1];
Rv2 = [cosd(-ap) 0 -sind(-ap); 0 1 0; sind(-ap) 0 cosd(-ap)];
V = Rv2*Rv1*V

n = [cosd(foraz)*sind(inc); sind(foraz)*sind(inc); cosd(inc)]

V*n
% end

