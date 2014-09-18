function [ inc ] = rayp2inc( rayp,V,r )
% [ inc ] = rayp2inc( rayp,v,r )
% Convert a ray parameter in seconds per degree to an incidence angle in 
% degrees, using the velocity, V (km/s) and a radius, r (km)
% IF no depth is given, assume this is a flat Earth and rayp is in s/km

if nargin < 3
    r = 6371;
end
z = 6371-r;

if nargin > 2
    inc = asind( (r2d(rayp) .* V)./(6371-z) );
else
    inc = asind( rayp.*V );
end

end

