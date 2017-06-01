function skm = rayp_sdeg2skm( sdeg,depth )

if nargin < 2 || isempty(depth)
    depth=0;
end

skm = sdeg * 360/(2*pi*(6371-depth));

end

