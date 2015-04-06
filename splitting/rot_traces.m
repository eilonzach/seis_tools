function [ datZYX_ ] = rot_traces( datZYX, theta )
% [ datZYX_ ] = rot_traces( datZYX, theta )
%   Function to rotate traces from right handed coordinates (ZYX) to new
%   right handed coordinates (ZY'X') where Z is positive downwards
% 
% INPUTS
% datZYX    - Nx3(2) matrix of seismogram traces in the order Z,Y,X (Y,X)
% theta     - clockwise angle through which to rotate traces
% 
% OUTPUTS
% datZYX_    - Nx3(2) matrix of seismogram traces in the order Z,Y',X' (Y,X)

R = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]; %Rotation matrix

if size(datZYX,1)< size(datZYX,2), datZYX = datZYX'; end % flip if in rows
datZYX_ = zeros(size(datZYX));

if size(datZYX,2)==2 % only horiz components
datZYX_ = datZYX*R';
elseif size(datZYX,2)==3 % 3-component
datZYX_(:,1) = datZYX(:,1);
datZYX_(:,2:3) = datZYX(:,2:3)*R';
end

end

