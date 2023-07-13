function [stainfo] = stainfo_extractfromIRIS_resp(stations_IRIS,stainfo)
%R[stainfo] = stainfo_extractfromIRIS_resp(stations_IRIS,stainfo)
%   This function is designed to extract the response information from a
%   station_INFO object obtained from an irisFetch.Stations('RESPONSE'....)
%   request. 

%% first extract basic station information if not already in existence
if ~exist('stainfo','var') || isempty(stainfo)

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

end % on whether need to build stainfo

%% Now move on to response information
if isempty(stations_IRIS(1).Channels(1).Response.Stage)
    warning('Response field seems to be empty')
    warning('Check irisFetch response at "RESPONSE" level')
    return
end

% start populating response info
stainfo.chanresps = cell(size(stainfo.chans));
% loop through stations and channels
% for each, concatenate poles, zeros, gain, input units, output units

for iso = 1:stainfo.nstas  % this is station in db
    ista = find(strcmp({stations_IRIS.StationCode},stainfo.stas(iso)) &...
                strcmp({stations_IRIS.NetworkCode},stainfo.nwk(iso)));
    if length(ista) > 1
%         keyboard
    end
    % make structure to house response
    resp_sta = cell(1,stainfo.nchans(iso)); % one cell element per channel. 
    % responses will be stored in structures, with one layer of structure
    % per response file per channel
    for ico = 1:stainfo.nchans(iso) % this is channel in db
        chan = stainfo.chans{iso,ico};
        % prep response structure for this channel
        resp_cha = cell2struct({1,[],[],nan,nan,'','','',''},{'Gain','Poles','Zeros','StartDate_ser','EndDate_ser','StartDate_str','EndDate_str','inUnits','outUnits'},2);
        % okay, now we know which station and which channel, time to go
        % through IRIS-obtained object and get all the responses
        jresp = 1;
        for isi = ista % loop through station numbers
            for ici = 1:length(stations_IRIS(isi).Channels) % loop through stations
                if ~strcmp(stations_IRIS(isi).Channels(ici).ChannelCode,chan), continue; end
                if isempty(stations_IRIS(isi).Channels(ici).Response.Stage), 
                    %No response in here
                    continue; 
                end
                % if you got here, it's a matching sta and channel. 
                resp_cha(jresp).StartDate_str = stations_IRIS(isi).Channels(ici).StartDate;
                resp_cha(jresp).EndDate_str   = stations_IRIS(isi).Channels(ici).EndDate;
                if isempty(resp_cha(jresp).EndDate_str)
                    resp_cha(jresp).EndDate_str = '2599-12-31 23:59:59.999';
                end
                resp_cha(jresp).StartDate_ser = datenum(resp_cha(jresp).StartDate_str);
                resp_cha(jresp).EndDate_ser   = datenum(resp_cha(jresp).EndDate_str);
                if isempty(stations_IRIS(isi).Channels(ici).CalibrationUnits)
                    resp_cha(jresp).inUnits = 'UNKNOWN';
                else
                    resp_cha(jresp).inUnits = stations_IRIS(isi).Channels(ici).CalibrationUnits.Name;
                end
                resp_cha(jresp).outUnits = stations_IRIS(isi).Channels(ici).SensitivityUnits;
                resp_cha(jresp).Gain = 1; % default
                resp_cha(jresp).Poles = [];
                resp_cha(jresp).Zeros = [];
                % the following is a response for the right station and
                % channel, now just have to work out date and values
                R = stations_IRIS(isi).Channels(ici).Response;
                % loop through response stages
                for istg = 1:length(R.Stage)
                    % multiply through gain
                    resp_cha(jresp).Gain = resp_cha(jresp).Gain*R.Stage(istg).StageGain.Value;
                    % if FIR, skip
                    if ~isempty(R.Stage(istg).FIR)
                        continue; 
                    end
                    % otherwise, must be poles and zeros
                    % concat on poles, zeros
                    if ~isempty(R.Stage(istg).PolesZeros)
                        resp_cha(jresp).Zeros = sort([resp_cha(jresp).Zeros, R.Stage(istg).PolesZeros.Zero]);
                        resp_cha(jresp).Poles = sort([resp_cha(jresp).Poles, R.Stage(istg).PolesZeros.Pole]);
                    end
                end
                % check gain is within 1% of "Sensitivity Value"
                if abs(1 - resp_cha(jresp).Gain./R.InstrumentSensitivity.Value) >0.01 
                    warning('Gain seems to be wrong')
                    fprintf('integrated gain for %s-%s %s between %s and %s is \n %.4f\n',...
                        stainfo.stas{iso},stainfo.nwk{iso},stainfo.chans{iso,ico},...
                        resp_cha(jresp).StartDate_str,resp_cha(jresp).EndDate_str,resp_cha(jresp).Gain)
                    fprintf('versus Instrument Sensitivity value of \n %.4f\n',...
                        R.InstrumentSensitivity.Value)
%                     keyboard
                else 
                    fprintf('DONE %s-%s %s between %s and %s \n',...
                        stainfo.stas{iso},stainfo.nwk{iso},stainfo.chans{iso,ico},...
                        resp_cha(jresp).StartDate_str,resp_cha(jresp).EndDate_str)
                    
                end
                % done with this response, move to next
                jresp = jresp + 1;
            end % loop on ici irisFetch chan
        end % loop on isi irisFetch station 
        % done with this db channel
        resp_sta{ico} = resp_cha;
    end % loop on ico chans for this db station
    stainfo.chanresps(iso,1:stainfo.nchans(iso)) = resp_sta;
end % loop on iso stas in db

    




end

