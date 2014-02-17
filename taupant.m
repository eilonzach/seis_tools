function [ tt ] = taupant( varargin )
%function to do similar things to tauptime (etc.) but using the antelope
%phase arrival time predictions and slownesses
%
% INPUTS:    taupant(...)
%             tt=taupant(...)
%             tt=taupant(...,'z|h|dep|edep|depth',depth,...)
%             tt=taupant(...,'p|ph|phase|iphase',phases,...)
%             tt=taupant(...,'d|deg|gcarc|degrees|delta',gcarc,...)
%             tt=taupant(...,'s|st|sta|station',[slat slon],...)
%             tt=taupant(...,'e|ev|evt|event',[elat elon],...)
% 
% default edep is 0 km
%  
% OUTPUTS:   a structure array with the following fields:
%               tt(index).time         - travel time (sec)
%                        .depth        - depth of earthquake (km)
%                        .distance     - true distance (deg)
%                        .rayparameter - ray parameter (sec/deg)
%                        .phase        - seismic phase name
%                        .station      - stn location if known ([lat lon])
%                        .event        - event location if known ([lat lon])
%                        .seaz         - station-event azimuth/backaz (deg)
%                        .esaz         - event-station azimuth/foraz  (deg)
%                        
%       Each phase has its own indice in the struct array TT.  Use TT(index)
%       to access individual phase information.


%% Default values
edep = 0;

%% Varargin
for i = 1:2:nargin
if(isempty(varargin{i})); continue; end
    switch lower(varargin{i})
        case {'s' 'st' 'sta' 'station'}
            slat = varargin{i+1}(1);
            slon = varargin{i+1}(2);
        case {'e' 'ev' 'evt' 'event'}
            elat = varargin{i+1}(1);
            elon = varargin{i+1}(2);
        case {'z' 'h' 'dep' 'edep' 'depth'}
            edep = varargin{i+1};
        case {'d' 'deg' 'gcarc' 'degrees' 'delta'}
            gcarc = varargin{i+1};
        case {'p','ph' 'phase' 'iphase' 'phases'}
            phase = varargin{i+1};
    end
end

%% Calculate needed values
if ~exist('gcarc','var')
    if exist('elat','var') && exist('slat','var')
        [gcarc,esaz] = distance(elat,elon,slat,slon);
        [kmlen,seaz] = distance_km(slat,slon,elat,elon);
        kmpd = kmlen/gcarc;
    else
        error('Need either distance or station and event to be specified')
    end
else
        kmpd = 111.1330;
end
N = length(phase);


%% times and slownesses
[artimes,phasenames] = arrtimes(gcarc,edep);
[slows,~] = arr_slowness(gcarc,edep); % in seconds/km
slows = slows*kmpd;

if N > 2
    pp = find(strncmp(phasenames,phase,N));
else
    pp = find(strcmp(phasenames,phase));
end

if isempty(pp)
    fprintf('NO matching phases found at this distance\n');
    tt = [];
    return
end

for ii = 1:length(pp)
    ip = pp(ii);
tt(ii) = struct('time',artimes(ip),'depth',edep,'distance',gcarc,...
                    'rayparameter',slows(ip),'phase',phasenames(ip),...
                    'station','','event',nan,'seaz',nan,'esaz',nan);

if exist('elat','var') && exist('slat','var')
    tt(ii).station = [elat,elon];
    tt(ii).event = [slat,slon];
    tt(ii).seaz = seaz;
    tt(ii).esaz = esaz;
end
end

end

