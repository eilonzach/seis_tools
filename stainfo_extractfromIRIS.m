function [stainfo] = stainfo_extractfromIRIS(stations_IRIS)
% [stainfo] = stainfo_extractfromIRIS(stations_IRIS)
%   This function is designed to extract the station information from a
%   station_INFO object obtained from an irisFetch.Stations('CHANNEL'....)
%   request. 

% find and fill in empty enddates
for is = 1:length(stations_IRIS)
    if isempty(stations_IRIS(is).EndDate)
        stations_IRIS(is).EndDate = '2599-12-31 23:59:59.999';
    end
end

stainfo = struct('stas',{{stations_IRIS.StationCode}'},...
                 'nwk',{{stations_IRIS.NetworkCode}'},...
                 'slats',[stations_IRIS.Latitude]',...
                 'slons',[stations_IRIS.Longitude]',...
                 'selevs',[stations_IRIS.Elevation]',...
                 'ondate',datenum({stations_IRIS.StartDate}'),...
                 'offdate',datenum({stations_IRIS.EndDate}'),...
                 'ondate_str',{{stations_IRIS.StartDate}'},...
                 'offdate_str',{{stations_IRIS.EndDate}'},...
                 'nstas',length(stations_IRIS));  
             

% parse channels             
chans = cell(stainfo.nstas,3);
chandips = nan(stainfo.nstas,3);
chanazs = nan(stainfo.nstas,3);
nchans = zeros(stainfo.nstas,1);
for is = 1:stainfo.nstas
    nchan = length(stations_IRIS(is).Channels);
    tempchans = cell(1,nchan);
    tempdips = nan(1,nchan);
    tempazs = nan(1,nchan);
    for ic = 1:nchan
        tempchans(ic) = {stations_IRIS(is).Channels(ic).ChannelCode};
        tempdips(ic) = stations_IRIS(is).Channels(ic).Dip;
        tempazs(ic) = stations_IRIS(is).Channels(ic).Azimuth;
    end
    [stachans,indch] = unique(tempchans);
    nchans(is) = length(stachans);
    chans(is,1:nchans(is)) = stachans;
    chandips(is,1:nchans(is)) = tempdips(indch);
    chanazs(is,1:nchans(is)) = tempazs(indch);
end
chandips(cellfun('isempty',chans)) = nan;
chanazs(cellfun('isempty',chans)) = nan;
stainfo.nchans = nchans;
stainfo.chans = chans;
stainfo.chandips = chandips;
stainfo.chanazs = chanazs;

[stainfo] = stainfo_unique(stainfo);






end

