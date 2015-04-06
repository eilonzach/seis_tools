function [ datZRT ] = zne2zrt( datZNE, baz )
% [ datZRT ] = zne2zrt( datZNE, baz )
%   Function to rotate traces from geographic coordinates (ZNE) to ray
%   coordinates (ZRT)
% 
% INPUTS
% datZNE    - a Nx3(2) matrix of seismogram traces in the order Z,N,E (N,E)
% baz       - backazimuth of arrival (in degrees)
% 
% OUTPUTS
% datZRT    - a Nx3(2) matrix of seismogram traces in the order Z,R,T (R,T)

foraz = mod(baz+180,360);
R = [cosd(foraz) sind(foraz); -sind(foraz) cosd(foraz)]; %Rotation matrix

if size(datZNE,1)< size(datZNE,2), datZNE = datZNE'; end % flip if in rows
datZRT = zeros(size(datZNE));

if size(datZNE,2)==2 % only horiz components
datZRT = datZNE*R';
elseif size(datZNE,2)==3 % 3-component
datZRT(:,1) = datZNE(:,1);
datZRT(:,2:3) = datZNE(:,2:3)*R';
end

end

