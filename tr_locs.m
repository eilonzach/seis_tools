function [locs,nlocs] = tr_locs(tr)
% locs = tr_locs(tr) 
%
% Function to find how many locations there are per trace (i.e. different
% seismometers/depths) assuming "tr" is a irisFetch-style trace object.

% assume one location
locs = ones(length(tr),1);
nlocs = 1;

% lalo
if length(unique([tr.latitude] + 1i*[tr.longitude])) > 1
    fprintf('Multiple physical locations for this station\n')
    phys_locs = unique([tr.latitude] + 1i*[tr.longitude]);
    nplocs = length(phys_locs);
    iloc = 1; 
    for ipl = 1:nplocs
        ind = [tr.latitude] + 1i*[tr.longitude] == phys_locs(ipl);
        if ~any(ind), continue; end
        locs(ind) = iloc;
        iloc = iloc+1;
    end
    nlocs = length(unique(locs));
end

% depth
if length(unique([tr.depth])) > 1
    fprintf('Multiple depths for this station\n')
    depth_locs = unique([tr.depth]);
    ndlocs = length(depth_locs);
    iloc = 1; loctemp = locs;
    for il = 1:nlocs
    for idl = 1:ndlocs
        ind = [tr.depth]' == depth_locs(idl) & locs(:) == il;
        if ~any(ind), continue; end
        loctemp(ind) = iloc;
        iloc = iloc+1;
    end
    end
    locs = loctemp;
    nlocs = length(unique(locs));    
end

% location
if length(unique([tr.location])) > 1
    fprintf('Multiple locations for this station\n')
    loc_locs = unique({tr.location});
    nllocs = length(loc_locs);
    iloc = 1; loctemp = locs;
    for il = 1:nlocs
    for ill = 1:nllocs
        ind = strcmp({tr.location}',loc_locs(ill)) & locs(:) == il;
        if ~any(ind), continue; end
        loctemp(ind) = iloc;
        iloc = iloc+1;
    end
    end
    locs = loctemp;
    nlocs = length(unique(locs));    
end

    

end