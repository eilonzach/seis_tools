function [M0, str, dip, rak] = momenttensor_invert(M,coordsystem)
%[M0, str, dip, rak] = momenttensor_invert(M,coordsystem)
%   This function inverts a 3x3 moment tensor for the two nodal planes,
%   each described by their strike, dip, and rake. 
% 
% INPUTS 
%   M - 3x3 moment tensor
%   coordsystem - "rtp" [default] for the mode coordinate system 
%                   (radial, -colat, long) i.e. (Up, South, East)
%               - "nwz" for the geographic coord system (North, West, Up)
%               - "nw_z" for the geographic coord system (North,West,Down)

%  OUTPUTS
%   M0 - Scalar moment
%   str - 2x1 vector of strikes (degrees) for planes 1 and 2 respectively
%   dip - 2x1 vector of dips (degrees) for planes 1 and 2 respectively
%   rak - 2x1 vector of rakes (degrees) for planes 1 and 2 respectively
% 
%    Mrtp = [Mrr Mrt Mrp       Mnwz = [Mnn Mnw Mnz
%            Mrt Mtt Mtp               Mnw Mww Mwz
%            Mrp Mtp Mpp];             Mnz Mwz Mzz];
% 
%       Z. Eilon	21 May 2013     

if nargin<2
    coordsystem = 'rtp';
end

if strcmp(coordsystem,'rtp')
R = [ 0    -1     0 % rotation matrix from Up,S,E to N,W,Down
      0     0    -1
     -1     0     0];  
elseif strcmp(coordsystem,'nwz')
R = [1,0,0;0,1,0;0,0,-1]; % rotation matrix from N,W,Up to N,W,Down  
elseif strcmp(coordsystem,'nw_z')
R = eye(3);
end

Mnw_z = R*M*R'; % moment tensor in (North,West,Down) coordinate system
[V,L] = eig(Mnw_z);
L = diag(L);
M0 = mean(abs(L([1,3])));

s = sqrt(0.5)*(V(:,3) + V(:,1));
n = sqrt(0.5)*(V(:,3) - V(:,1));
b = cross(n,s); % null vector
%% First nodal plane
dip(1) = acosd(n(3));
rak(1) = r2d(atan2(s(3),b(3)));
str(1) = mod(r2d(atan2(n(1),n(2))),360);
if dip(1)>90
dip(1) = 180-dip(1);
rak(1) = 360-rak(1);
str(1) = 180+str(1);
end
    
%% Second plane 
% same as first, but switch s and n
dip(2) = acosd(s(3));
rak(2) = r2d(atan2(n(3),-b(3)));
str(2) = mod(r2d(atan2(s(1),s(2))),360);
if dip(2)>90
dip(2) = 180-dip(2);
rak(2) = 360-rak(2);
str(2) = 180+str(2);
end
str = mod(str,360);
rak = mod(rak+180,360)-180;

fprintf('Scalar Moment = %.2e\n',M0)
fprintf('Fault plane:  strike=%.0f \t dip=%.0f \t slip=%.0f\n',str(1),dip(1),rak(1))
fprintf('Fault plane:  strike=%.0f \t dip=%.0f \t slip=%.0f\n',str(2),dip(2),rak(2))

end


