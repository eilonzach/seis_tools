function [ inc ] = rayp2inc( rayp,V,z )
% [ inc ] = rayp2inc( rayp,v,z )
% Convert a ray parameter in seconds per degre to an incidence angle, using
% the velocity, V (km/s) and a depth, z (km)
% IF no depth is given, assume this is a flat Earth and rayp is in s/km

if nargin>2
    inc = asind( (r2d(rayp) .* V)./(6371-z) );
end
if nargin < 3
    inc = asind( rayp.*V );
end

end

