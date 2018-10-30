function [stainfo] = stainfo_unique(stainfo)
%[stainfo] = stainfo_unique(stainfo)
%   
% function to parse the stainfo structure to remove any repeat stations

% get unique station and network names
unqstas = unique([stainfo.stas]); nustas = length(unqstas);
unqnwk = unique([stainfo.nwk]);   nunwk = length(unqnwk);

kill = false(stainfo.nstas,1);
for in = 1:nunwk
for is = 1:nustas
    inds = find(strcmp([stainfo.stas],unqstas(is)) & strcmp([stainfo.nwk],unqnwk(in)));
    if length(inds)>1
        kill(inds(1:end-1)) = true;
        [stainfo.ondate(inds(end)),ibeg]  = min(stainfo.ondate(inds));
        [stainfo.offdate(inds(end)),iend] = max(stainfo.offdate(inds));
        stainfo.ondate_str(inds(end))  = stainfo.ondate_str(inds(ibeg));
        stainfo.offdate_str(inds(end)) = stainfo.offdate_str(inds(iend));
    end
end
end

fns = fieldnames(stainfo);
for iff = 1:length(fns)
    if length(stainfo.(fns{iff}))>1
        stainfo.(fns{iff})(kill) = [];
    end
end
stainfo.nstas = length(stainfo.stas);
    

end

