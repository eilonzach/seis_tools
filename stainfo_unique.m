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
        % keep last one of multiple, kill the rest
        kill(inds(1:end-1)) = true;
        % set beginning and end dates to biggest span
        [stainfo.ondate(inds(end)),ibeg]  = min(stainfo.ondate(inds));
        [stainfo.offdate(inds(end)),iend] = max(stainfo.offdate(inds));
        stainfo.ondate_str(inds(end))  = stainfo.ondate_str(inds(ibeg));
        stainfo.offdate_str(inds(end)) = stainfo.offdate_str(inds(iend));
        % if channels, concat all, 
        stachans    = stainfo.chans(inds,:);
        stachanazs  = stainfo.chanazs(inds,:);
        stachandips = stainfo.chandips(inds,:);
        % find unique chans
        istachans = find(~cellfun(@isempty,stachans));
        [~,iustachans] = unique(stachans(istachans));
        % unique names, azs, dips
        stachans    = stachans(istachans(iustachans));
        stachanazs  = stachanazs(istachans(iustachans));
        stachandips = stachandips(istachans(iustachans));
        % reinsert
        stainfo.nchans(inds(end)) = length(stachans);
        stainfo.chans(inds(end),1:stainfo.nchans(inds(end))) = stachans;
        stainfo.chanazs(inds(end),1:stainfo.nchans(inds(end))) = stachanazs;
        stainfo.chandips(inds(end),1:stainfo.nchans(inds(end))) = stachandips;
    end
end
end

fns = fieldnames(stainfo);
for iff = 1:length(fns)
    if length(stainfo.(fns{iff}))>1
        stainfo.(fns{iff})(kill,:) = [];
    end
end
stainfo.nstas = length(stainfo.stas);
    

end

