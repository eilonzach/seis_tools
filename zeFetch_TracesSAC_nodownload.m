function traces = zeFetch_TracesSAC_nodownload(nwk,sta,loc,chan,time_start,time_end,includePZ,predownloaddir)
% traces = zeFetch_TracesSAC_nodownload(nwk,sta,loc,chan,time_start,time_end,includePZ,predownloaddir)
%
% Note, this function has input and output identical to zeFetch_Traces.m,
% but instead of downloading and reading mseed files, it downloads and
% reads SAC files from the EarthScope DMC, which might be a more stable
% format.
%
% INPUTS
%  nwk  = two-character code
%  sta  = four(or five) character code
%  loc = '*' % two char code
%  chan = '?H?'; % three char code
%  time_start = '2010-02-27T06:30:00.000' % string
%  time_end = '2010-02-27T07:30:00.000' % string
%  includePZ = true; % logical
%
% Z. Eilon 03/2025

arguments
    nwk
    sta
    loc = '*';
    chan = '?H?';
    time_start = '2010-02-27T06:30:00.000'
    time_end = '2010-02-27T07:30:00.000'
    includePZ = true;
    predownloaddir = '.';
end

spd = 24*3600;

addpath('~/Dropbox/MATLAB/lib/sacutils/')

% put times in right format
time_start = regexprep(time_start,' ','T');
time_end = regexprep(time_end,' ','T');
time_start_serial = convertTo(datetime(time_start),'datenum');
time_end_serial = convertTo(datetime(time_end),'datenum');

%% Do the data request
% build URL for data request
data_file_searchstring = sprintf('%s.%s.*.SAC',nwk,sta);

sacfiles = dir(fullfile(predownloaddir,data_file_searchstring));

if isempty(sacfiles)
    error('%s\nNo data seems to be available - empty file returned',data_file_searchstring)
    return
end

%% Do the response request
% build URL for response request
resp_url_preamble = 'http://service.iris.edu/irisws/sacpz/1/query?';
resp_url = [resp_url_preamble,...
       sprintf('net=%s&',nwk),...
       sprintf('sta=%s&',sta),...
       sprintf('loc=%s&',loc),...
       sprintf('cha=%s&',chan),...
       sprintf('start=%s&',time_start),...
       sprintf('end=%s',time_end),...
        ];
% e.g. 'http://service.iris.edu/fdsnws/dataselect/1/query?net=IU&sta=ANMO&loc=00&cha=BHZ&start=2010-02-27T06:30:00.000&end=2010-02-27T10:30:00.000'

% get responses as pz
pz=readsacpz_rdseed(urlread(resp_url),true);
% small tweak to get start and end times in serial date
for ipz = 1:length(pz.path)
    pz.b_str{ipz,1} = sprintf('%4.0f-%03.0f %02.0f:%02.0f:%02.0f',pz.b(ipz,:));
    pz.b_dt(ipz,1) = datetime(pz.b_str{ipz},'Format','uuuu-DDD HH:mm:ss');
    pz.e_str{ipz,1} = sprintf('%4.0f-%03.0f %02.0f:%02.0f:%02.0f',pz.e(ipz,:));
    pz.e_dt(ipz,1) = datetime(pz.e_str{ipz},'Format','uuuu-DDD HH:mm:ss');
end

%% now read the data
addpath('~/Dropbox/MATLAB/lib/sacutils/');

% parse filename
sacfileinfo = parse_sac_filename({sacfiles.name});

%% prepare traces object
traces0 = struct('network',[],... y
                'station',[],... y
                'location','',... y
                'channel',[],... y
                'quality',[],... y
                'latitude',[],...
                'longitude',[],...
                'elevation',[],...
                'depth',[],...
                'azimuth',[],...
                'dip',[],...
                'sensitivity',[],...
                'sensitivityFrequency',[],...
                'instrument',[],...
                'sensitivityUnits',[],...
                'data',[],... y
                'sampleCount',[],... y
                'sampleRate',[],... y
                'startTime',[],... y
                'endTime',[],... y
                'sacpz',struct('units',[],'constant',[],'poles',[],'zeros',[])...
                );

%% determine how many we are stations dealing with
stas = unique({sacfileinfo.sta});
Nstas = length(stas);

%% now for each station, loop through channels. pull out data
% set up empty traces object, will fill.
traces = traces0;
itr = 0;

for is = 1:Nstas
    sta = stas{is};
    ista = strcmp(sta,{sacfileinfo.sta});
    chans = unique({sacfileinfo(ista).chan});
    Nchans = length(chans);
    for ic = 1:Nchans        
        chan = chans{ic};
        ichan = ista & strcmp(chan,{sacfileinfo.chan});
        locs = unique({sacfileinfo(ichan).loc});
        Nlocs = length(locs);
        for il = 1:Nlocs
            loc = locs{il};
            iloc = ichan & strcmp(loc,{sacfileinfo.loc});

            sacfile = fullfile(sacfiles(iloc).folder, sacfiles(iloc).name);
            %% HERE ACTUALLY PULL OUT DATA FROM SAC FILE
            sinfo = readsac(sacfile); % get all metadata, together with data
            % [stt,svv] = readsac(sacfile); % get time and velocity vectors

            % % might have to patch due to odd mseed breakups
            % if length(b) > 1
            %     b = patch_gaps(b,time_start_serial,time_end_serial - 1/b(1).sampleRate/spd);
            % end
            % % might have to shave off last sample if an extra odd one
            % try
            %     b = shave_or_backfill_end(b,time_start_serial,time_end_serial);
            % catch
            %     try
            %         b = patchzeros(b,time_start_serial,time_end_serial);
            %         b = shave_or_backfill_end(b,time_start_serial,time_end_serial);
            %     catch
            %         keyboard
            %     end
            % end

            % now assign data into traces object
            itr = itr+1;
            traces(itr,1).network     = sinfo.KNETWK;
            traces(itr,1).station     = sta;
            traces(itr,1).location    = loc;
            traces(itr,1).channel     = chan;
            % traces(itr,1).quality     = b.dataquality;  % not sure if this is correct.
            traces(itr,1).sampleRate  = round(abs(1./sinfo.DELTA));
            traces(itr,1).sampleCount = sinfo.NPTS;
            traces(itr,1).data        = double(sinfo.DATA1(:));  % columnise, double-ize

            % get times exactly
            [NZMO,NZDAY] = calday(sinfo.NZYEAR,sinfo.NZJDAY);

            traces(itr,1).startTime   = datenum([sinfo.NZYEAR,NZMO,NZDAY,sinfo.NZHOUR,sinfo.NZMIN,sinfo.NZSEC+sinfo.NZMSEC/1000]); 
            tr_dur = seconds((traces(itr,1).sampleCount-1)./traces(itr,1).sampleRate);
            traces(itr,1).endTime     = traces(itr,1).startTime + days(tr_dur);
            % let's put in datetime format, too
            traces(itr,1).startTime_dt = datetime([sinfo.NZYEAR,NZMO,NZDAY,sinfo.NZHOUR,sinfo.NZMIN,sinfo.NZSEC+sinfo.NZMSEC/1000]);
            traces(itr,1).endTime_dt   = traces(itr,1).startTime_dt + tr_dur;
    
            % now kill any data after time_end or before time_start
            tt = traces(itr,1).startTime_dt:seconds(1./traces(itr,1).sampleRate):traces(itr,1).endTime_dt;
            oktime = tt >= datetime(time_start) & tt <= datetime(time_end);
            ttok = tt(oktime);
            % SHOULD ADD FUNCTION TO  zero-pad if needed

            traces(itr,1).startTime_dt = ttok(1);
            traces(itr,1).endTime_dt = ttok(end);
            traces(itr,1).data = traces(itr,1).data(oktime);
            traces(itr,1).sampleCount = sum(oktime);
            traces(itr,1).startTime = datenum(traces(itr,1).startTime_dt);
            traces(itr,1).endTime = datenum(traces(itr,1).endTime_dt);


            % now get the rest of the information from the pz object
            ipz = strcmp(sta,pz.kstnm) & strcmp(chan,pz.kcmpnm);
            if ipz==0, continue; end % NO RESP INFORMATION FOR THIS CHANNEL
            if sum(ipz)>1
                ipz = ipz & (pz.e_dt > traces(itr).endTime_dt) & (pz.b_dt < traces(itr).startTime_dt); % limit by time of resp file
            end
            if sum(ipz)>1
                % choose the most recent 
                ipzs = find(ipz);
                [~,jj] = max(pz.b_dt(ipzs));
                ipz = ipzs(jj);
                % keyboard % really should be unique response at this point
            end
    
            traces(itr,1).latitude    = pz.stla(ipz);            
            traces(itr,1).longitude   = pz.stlo(ipz);
            % check equality to the SAC info
            dd = (traces(itr,1).latitude - sinfo.STLA).^2 + ((traces(itr,1).longitude - sinfo.STLO)*cosd(traces(itr,1).latitude)).^2;
            if dd > (0.01/111)^2
                error('PZ and the SAC file seem to have location differing by >10m');
            end

            traces(itr,1).elevation   = pz.stel(ipz);
            traces(itr,1).depth       = pz.stdp(ipz);
            traces(itr,1).azimuth     = pz.cmpaz(ipz);
            traces(itr,1).dip         = 90-pz.cmpinc(ipz); % go from inclination to dip (former is measured from vertical, latter from horizontal.
            traces(itr,1).input       = pz.input{ipz}; % not sure if this is correct.
            traces(itr,1).sensitivity = pz.sensitivity(ipz);
            traces(itr,1).instrument  = pz.insttype{ipz};
            traces(itr,1).sensitivityUnits     = pz.sensitivityunits{ipz};
            traces(itr,1).sensitivityFrequency = [];
            if includePZ
            traces(itr,1).sacpz.units     = pz.instgainunits{ipz};
            traces(itr,1).sacpz.constant  = pz.k(ipz);
            traces(itr,1).sacpz.poles     = pz.p{ipz};
            traces(itr,1).sacpz.zeros     = pz.z{ipz};
            traces(itr,1).sacpz.A0        = pz.a0(ipz);
            % attempt to parse if extra zero
                ezero = false;
                if strcmp(pz.input{ipz},'M') && ...
                        ... next line tests if by 1Hz there is one more zero than pole (should be same if flat here)
                        (sum(abs(pz.z{ipz}) < 1) - sum(abs(pz.p{ipz}) < 1))==1 
                    ezero = 1;
                end
            end
            traces(itr,1).sacpz.extrazerotoV = ezero;

            % DELETE unneeded file
        end % on locations
    end % on channels
end % on stations

end



%% subfunctions
function a_patched = patch_gaps(a,time_start,time_end)

if nargin <2
    time_start = a(1).matlabTimeVector(1);
    time_end   = a(end).matlabTimeVector(end);
end

spd = 24*3600;
dt = abs(1./a(1).sampleRate);
% first try to do very simplistically - just concatenate to end
a_patched = a(1);
for ii = 2:length(a)
    a_patched.data = [a_patched.data;a(ii).data];
    a_patched.matlabTimeVector = [a_patched.matlabTimeVector;a(ii).matlabTimeVector];
    a_patched.sampleCount = a_patched.sampleCount+a(ii).sampleCount;
end
% test if this is adequate
if round(max(diff(a_patched.matlabTimeVector))*spd,3) == dt
    return % yes, concatenation solved it. end
end
% if we get to this point, some samples missing. Must interpolate. 
% Might be big gaps. Might be some reversed times?!
% Might also be some discontinuities in waveform

% find indices at which time vector seems to jump. This might separate fine
% but improperly stitched sections of data. 
dt_bad = find(round_level(abs(diff(a_patched.matlabTimeVector - a_patched.matlabTimeVector(1))*spd),0.1*dt) > 1.01*dt);
Ntbad = length(dt_bad);
% note ^ are indices at the beginning of the pairs that are a time jump
% from plotting, seems like there are some discontinuities in the data
% around some/all of these time jumps. 
% find associated data discontinuities. 
dd_reljump = zeros(Ntbad,1);
% loop through and see if the data jump at the dt_gap is greater than the
% local average
for ibad = 1:Ntbad
    % find indices of 100 samples on either side of the gap. We'll use the
    % local diffs in this region
    j11 = max([1,dt_bad(ibad)-100]); % start of left 100 samples
    j12 = dt_bad(ibad); % end of left 100 samples
    j21 = dt_bad(ibad)+1; % start of right 100 samples
    j22 = min([a_patched.sampleCount,dt_bad(ibad)+100]); % end of right 100 samples
    if ibad < Ntbad
        j22 = min([j22,dt_bad(ibad+1)-100]); % in case we run into next bad dt to the right.
    end
    localmeandd = mean(abs([diff(a_patched.data(j11:j12));diff(a_patched.data(j21:j22))]));
    localstddd = std(abs([diff(a_patched.data(j11:j12));diff(a_patched.data(j21:j22))]));
    baddd = abs(diff(a_patched.data(j12:j21)));
    % this is the jump at the dt discontinuity, relative to the local mean,
    % expressed in multiples of standard deviations of the local value.
    % Anything above 2 is suspicious, 3 is v suspicious, 5 is bonkers. 
    dd_reljump(ibad) = (baddd - localmeandd)/localstddd; 
end

% now, consider patching the waveform by shifting segments
if any(dd_reljump>4)
    warning('Shifted data at anomalous jumps in trace associated with time errors')
    a_patched.warnflag = 'Shifted data at anomalous jumps in trace associated with time errors';
    newdat = a_patched.data;
        % while testing
        % figure(55); clf
        % stackedplot([1:a_patched.sampleCount]',[(a_patched.matlabTimeVector-a_patched.matlabTimeVector(1))*spd-[0:a_patched.sampleCount-1]'*dt,newdat]);
        % end while testin
    for ibad = 1:Ntbad
        j11 = max([1,dt_bad(ibad)-100]); % start of left 100 samples
        j12 = dt_bad(ibad); % end of left 100 samples
        j21 = dt_bad(ibad)+1; % start of right 100 samples
        j22 = min([a_patched.sampleCount,dt_bad(ibad)+100]); % end of right 100 samples
        if j21==j22 % there is a time gap at the very end...
            newdat(j21) = newdat(j12); continue;
        end
        localmeandd = mean(abs([diff(newdat(j11:j12));diff(newdat(j21:j22))]));
        localstddd = std(abs([diff(newdat(j11:j12));diff(newdat(j21:j22))]));
        baddd = abs(diff(newdat(j12:j21)));
        reljump = (baddd - localmeandd)/localstddd;
        proj21_from1 = 2*newdat(j12) - newdat(j12-1);
        proj12_from2 = 2*newdat(j21) - newdat(j21+1);
        dc_forward = newdat(j21) - proj21_from1;
        dc_backward = proj12_from2 - newdat(j12);
        dc = 0.5*(dc_forward + dc_backward);
        % this is the shift. Now move the next segment by this amount
        nseg1 = j21;
        if ibad < Ntbad
            nseg2 = dt_bad(ibad+1);
        else
            nseg2 = a_patched.sampleCount;
        end
        newdat(nseg1:nseg2) = newdat(nseg1:nseg2) - dc;
    end
        % while testing
        % figure(56); clf
        % stackedplot([1:a_patched.sampleCount]',[(a_patched.matlabTimeVector-a_patched.matlabTimeVector(1))*spd-[0:a_patched.sampleCount-1]'*dt,newdat]);
        % end while testin
    a_patched.data = newdat;
end


% check for linear time
if any(diff(a_patched.matlabTimeVector)<0)
    % argh, some times are reversed?
    warning('Nonsensical times; reset time axis assuming times are just plain wrong')
    a_patched.warnflag2 = 'Nonsensical times; reset time axis assuming times are just plain wrong';
    a_patched.matlabTimeVector = time_start + [0:a_patched.sampleCount-1].*dt/spd;
end

a_patched2 = a(1);
a_patched2.matlabTimeVector = [time_start:dt/spd:time_end]';
a_patched2.sampleCount = length(a_patched2.matlabTimeVector);
a_patched2.data = interp1(a_patched.matlabTimeVector,a_patched.data,a_patched2.matlabTimeVector,'linear','extrap');
a_patched2.data = a_patched2.data(:); % columnise
a_patched = a_patched2;

end


function a_shaved = shave_or_backfill_end(a,time_start_serial,time_end_serial)
spd = 24*3600;

nsamp_expected = round(spd*(time_end_serial - time_start_serial).* a(1).sampleRate);
if a.sampleCount == double(nsamp_expected)
    a_shaved = a; 
    return
else
    a.matlabTimeVector = round(a.matlabTimeVector,ceil(log10(spd*abs(a.sampleRate)))+3);
    if time_end_serial <= a.matlabTimeVector(end)
        sampkill = a.matlabTimeVector>=time_end_serial;
        a_shaved = a;
        a_shaved.data(sampkill) = [];
        a_shaved.matlabTimeVector(sampkill) = [];
        a_shaved.sampleCount = length(a_shaved.data);
        a_shaved.numberOfSamples = a_shaved.sampleCount;
    else
        warning('possibly not enough samples')
        a_shaved = a;
        a_shaved.data(a.sampleCount+1:nsamp_expected) = nan;
        nsamp_to_add = nsamp_expected-a.sampleCount;
        a_shaved.matlabTimeVector(a.sampleCount:nsamp_expected) = a_shaved.matlabTimeVector(a.sampleCount)+[0:nsamp_to_add]./a(1).sampleRate./spd;    
        a_shaved.sampleCount = length(a_shaved.data);
        a_shaved.numberOfSamples = a_shaved.sampleCount;        
    end
end


end



function a_patchzeros = patchzeros(a,time_start_serial,time_end_serial)
% to specifically deal with weird case where there seem to be a string of
% zeros at various places in the time series. 
spd = 24*3600;

samples_expected = round((time_end_serial-time_start_serial)*spd.*a.sampleRate);
samples_actual = a.sampleCount;

a_patchzeros = a;

% find the odd zeros
izero = find(a.data==0);
dizero = diff(izero);
nsegments = sum(dizero~=1)+1;
zero_segment_inds = cell(nsegments,1);
zero_segment_lens = nan(nsegments,1);
ind = 1;
for iseg = 1:nsegments
    while ind<=length(izero)
        zero_segment_inds{iseg} = [zero_segment_inds{iseg};izero(ind)]; 
        zero_segment_lens(iseg) = length(zero_segment_inds{iseg});
        ind = ind+1;
        if ind>length(izero), break; end
        if izero(ind)-izero(ind-1)~=1
            break
        end
    end
end
% check we got to the end
if sum(zero_segment_lens)~=length(izero)
    keyboard; 
end


% calculate jumps at segment ends as fractions of pre/post-ceeding 5 point
% jumps (basically looking for big discontinuities)
zero_segment_jumps = nan(nsegments,2);
for iseg = 1:nsegments
    lz = zero_segment_inds{iseg}(1);
    li = lz + [-6:-1]; 
    % account for silly case where gap is right at the beginning and no li exists
    li(li<1)=[]; 
    if numel(li>=2)
        lj_pre = mean(abs(diff(a.data(li))));
        lj = double(abs(a.data(lz) - a.data(li(end)))); 
        zero_segment_jumps(iseg,1) = lj/lj_pre;
    end
    rz = zero_segment_inds{iseg}(end);
    ri = rz + [1:6];
    % account for silly case where gap is right at the end and no/little ri exists
    ri(ri>a.numberOfSamples)=[]; 
    if numel(ri>=2)
        rj_post = mean(abs(diff(a.data(ri))));
        rj = double(abs(a.data(rz) - a.data(ri(1)))); 
        zero_segment_jumps(iseg,2) = rj/rj_post;
    end
end
zero_segment_jumps = nanmean(zero_segment_jumps,2);

% % now we have all the segments, throw out the segments of length 1. 
% segment_kill = zero_segment_lens==1;
% 
% if (samples_actual-samples_expected) == sum(zero_segment_lens(~segment_kill))
%     % grand, we have a match!! just kill everything in those segments
% else
%     % try killing segments with length of 2, esp if not first or last seg
%     segment_maybe_kill = zero_segment_lens==2;
%     segment_maybe_kill(1) = false; segment_maybe_kill(end) = false;
%     if (samples_actual-samples_expected) == sum(zero_segment_lens(~(segment_maybe_kill | segment_kill)))
%         % grand, we have a match!! just kill everything in those segments
%         segment_kill = segment_maybe_kill | segment_kill;
%     else
% 
%     % problem: the issue is more complicated than just some weird extra
%     % segments of zeros. 
%     keyboard
%     end
% end

% better way of doing it; entertain natural sequences of zeros up to length
% 3. But also don't allow anything with enormous jumps, even if small
% sequence

ok2kill = false;
for shortsegment = 1:3
    segment_kill = zero_segment_lens<=shortsegment & zero_segment_jumps<5;
    if (samples_actual-samples_expected) == sum(zero_segment_lens(~segment_kill))
        % grand, we have a match!! just kill everything in those segments
        ok2kill = true;
        break
    end
end
if ~ok2kill
    % problem: the issue is more complicated than just some weird extra
    % segments of zeros. 
    keyboard
    return
end



% do the killing
zero_segment_inds(segment_kill) = [];
izerokill = [];
for iseg = 1:length(zero_segment_inds)
    izerokill = [izerokill;zero_segment_inds{iseg}(:)];
end
a_patchzeros.data(izerokill) = [];
a_patchzeros.sampleCount = length(a_patchzeros.data);
a_patchzeros.numberOfSamples = length(a_patchzeros.data);

% if matlabTimeVector is blank....
if isempty(a_patchzeros.matlabTimeVector)
    t0 = datenum(a.dateTime.year,a.dateTime.month,a.dateTime.day,a.dateTime.hour,a.dateTime.minute,a.dateTime.second);
    tt = [0:(a_patchzeros.numberOfSamples-1)]'/a_patchzeros.sampleRate/spd + t0;
    a_patchzeros.matlabTimeVector = tt;
end




end