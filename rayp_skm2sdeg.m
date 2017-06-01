function sdeg = rayp_skm2sdeg( skm,depth )

if nargin < 2 || isempty(depth)
    depth=0;
end

sdeg = skm * 2*pi*(6371-depth)/360;

end

