function [goodevts]=gdevts_fn(startdate,enddate,MagLims,DepthLims,LatLims,LonLims,refpt,DeltaLims)
%[goodevts]=gdevts_fn(startdate,enddate,MagLims,DepthLims,LatLims,LonLims,refpt,DeltaLims)
%   like gdorids_fn, but instead of getting events that satisfy criteria
%   from a database, it gets it from the IRIS DMC
%
%% Set criteria for accepted events
% NB THESE VALUES ARE INCLUSIVE!! ( <= or >= )
% 
% TIME
% startdate = [start_y,start_m,start_d]; 
% enddate   = [end_y,  end_m,  end_d];   
% 
% MAGNITUDE
% MagLims   = [minMag maxMag];
% 
% LOCATION
% DepthLims = [minDepth maxDepth] (km)      [] for default: all
% LatLims   = [minLat maxLat]               [] for default: all
% LonLims   = [minLon maxLon]               [] for default: all
% 
% DISTANCE 
% refpt=[reflat,reflon]; % [lat,long] of reference point for distance constraint
% mindeg=85; %minimum angular distance (0?degrees?180)
% maxdeg=140; %maximum angular distance (0?degrees?180)


% Get times into right format
search_wstart = sprintf('%4.0f-%2.0f-%2.0f 00:00:00',startdate(1),startdate(2),startdate(3));
search_wend = sprintf('%4.0f-%2.0f-%2.0f 23:59:59',enddate(1),enddate(2),enddate(3));

if isempty(LatLims), LatLims = [-90 90]; end
if isempty(LonLims), LonLims = [-180 180]; end
if isempty(DepthLims), DepthLims = [0 1000]; end

%%  Get event information from IRIS database
evinfo = irisFetch.Events('boxcoordinates',[LatLims LonLims],...
		'startTime',search_wstart,'endTime',search_wend,...
		'minimumMagnitude',MagLims(1),'maximumMagnitude',MagLims(2),...
        'minimumDepth',DepthLims(1),'maximumDepth',DepthLims(2)');
    
if isempty(evinfo), error('No events match search criteria. Sorry...'); end
    
[evinfo.mag] = evinfo.PreferredMagnitudeValue; evinfo = rmfield(evinfo,'PreferredMagnitudeValue');
[evinfo.lat] = evinfo.PreferredLatitude; evinfo = rmfield(evinfo,'PreferredLatitude');
[evinfo.lon] = evinfo.PreferredLongitude; evinfo = rmfield(evinfo,'PreferredLongitude');
[evinfo.dep] = evinfo.PreferredDepth; evinfo = rmfield(evinfo,'PreferredDepth');
[evinfo.time] = evinfo.PreferredTime; evinfo = rmfield(evinfo,'PreferredTime');

% impose distance window
if exist('refpt','var')==1
    [dist,baz] = distance(refpt(1),refpt(2),[evinfo.lat],[evinfo.lon]);
    for ie = 1:length(evinfo)
        evinfo(ie).delta = dist(ie); 
        evinfo(ie).baz   = baz(ie); 
    end
    evinfo(dist > DeltaLims(2) | dist < DeltaLims(1)) = [];

else 
    for ie = 1:length(evinfo)
        evinfo(ie).delta = nan; 
        evinfo(ie).baz = nan; 
    end 
end


%% Tabulate
fprintf('ACTUALLY ONLY %u EVENTS MEET ALL CRITERIA\n\n',length(evinfo));
fprintf('Orid  Date       Time          Mag    Lat     Lon     Depth  Delta  Baz  Region\n\n')
for ie = 1:length(evinfo)
    fprintf('%4u  %s  %3.1f  %6.2f  %7.2f  %6.1f  %5.1f   %3.0f  %s\n',ie,...
        evinfo(ie).time,evinfo(ie).mag,...
        evinfo(ie).lat,evinfo(ie).lon,evinfo(ie).dep,...
        evinfo(ie).delta,evinfo(ie).baz,evinfo(ie).FlinnEngdahlRegionName)
end

goodevts = evinfo;

% plot


end
