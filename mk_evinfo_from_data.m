function [evinfo] = mk_evinfo_from_data(datadir,ofile)
% [evinfo] = mk_evinfo_from_data(datadir,ofile)
%   
% Function to get the event info from a data containing folders
% corresponding to earthquakes, with names "orid_yyyymmddHHMM"
% 
% Uses irisFetch to fill out other details. 
% 
% Option to specify a file to save the evinfo structure.

javaaddpath('/Users/zeilon/Dropbox/MATLAB/AttenBody/IRIS-WS-2.0.15.jar')



if nargin < 2
    ofile = [];
end

eqdirs = dir(datadir);

evinfo = struct('orids',0,...
            'norids',0,...
            'elats',0',...
            'elons',0',...
            'edeps',0',...
            'evmags',0',...
            'evtimes',0,...
            'evtimes_IRISstr',cell(0),...                 
            'datestamp',cell(0));                 


norids = 0;
for ifile = 1:length(eqdirs)
    if strcmp(eqdirs(ifile).name(1),'.')
        continue;
    end
    norids = norids+1;
    
    A = textscan(eqdirs(ifile).name,'%u_%4u%2u%2u%2u%2u');
    
    events_IRIS = irisFetch.Events('minimumMagnitude',6,...
                          'startTime',sprintf('%4u-%02u-%02u %02u:%02u:00',A{2},A{3},A{4},A{5},A{6}),...
                          'endTime',sprintf('%4u-%02u-%02u %02u:%02u:00',A{2},A{3},A{4},A{5},A{6}+1));
                      
	if length(events_IRIS)>1
        B = textscan(eqdirs(ifile+1).name,'%u_%4u%2u%2u%2u%2u');
        C = textscan(eqdirs(ifile-1).name,'%u_%4u%2u%2u%2u%2u');
        if events_IRIS(1).FlinnEngdahlRegionCode==events_IRIS(2).FlinnEngdahlRegionCode
            events_IRIS = events_IRIS(1);
        elseif isequal({A{2:end}},{B{2:end}}) % check if same minute as next one
                events_IRIS = events_IRIS(1); % keep only the first, not the next
        elseif isequal({A{2:end}},{C{2:end}}) % check if same minute as last one
                events_IRIS = events_IRIS(2); % keep only the second, not the last
        else
            events_IRIS = events_IRIS(mindex(-[events_IRIS.PreferredMagnitudeValue])); % just choose biggest
        end
    end

	evinfo(1).orids(norids,:) = double(A{1});
	evinfo.norids = norids;
    evinfo.elats(norids,:) = events_IRIS.PreferredLatitude;
    evinfo.elons(norids,:) = events_IRIS.PreferredLongitude;
    evinfo.edeps(norids,:) = events_IRIS.PreferredDepth;
    evinfo.evmags(norids,:) = events_IRIS.PreferredMagnitudeValue;
    evinfo.evtimes(norids,:) = datenum(events_IRIS.PreferredTime);
    evinfo.evtimes_IRISstr{norids,:} = events_IRIS.PreferredTime;           
    evinfo.datestamp{norids,:} = datestr(datenum(events_IRIS.PreferredTime),'yyyymmddHHMM');  
    
end

% sort into correct order so ie = orid
[~,sindex] = sort(evinfo.orids);  
evinfo.orids = evinfo.orids(sindex);
evinfo.elats = evinfo.elats(sindex);
evinfo.elons = evinfo.elons(sindex);
evinfo.edeps = evinfo.edeps(sindex);
evinfo.evmags = evinfo.evmags(sindex);
evinfo.evtimes = evinfo.evtimes(sindex);
evinfo.evtimes_IRISstr = evinfo.evtimes_IRISstr(sindex);
evinfo.datestamp = evinfo.datestamp(sindex);


if ~isempty(ofile)
    save(ofile,'evinfo');
end
