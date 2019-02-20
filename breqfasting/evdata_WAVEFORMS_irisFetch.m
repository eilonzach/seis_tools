function [datafile] = evdata_WAVEFORMS_irisFetch(station,network,request_details,respdir)
% [datafile] = evdata_WAVEFORMS_irisFetch(station,network,request_details,respdir)

if nargin < 1 || isempty(station) 
    station = 'J19A';
end
if nargin < 2 || isempty(network)
    network = 'TA';
end
if nargin < 3 || isempty(request_details)
    request_details = struct('phases',{{'P','S'}},'gclims',[65 75],'maglims',[5.7 7.2],...
                             'samprate',40,'datwind',[-100 100],'chans','HH?');
end

phases = request_details.phases;

datwind = request_details.datwind; % in sec around main arrival
samprate = request_details.samprate;
% tapertime = 5; % s at beginnning and end of window to taper over

if isfield(request_details,'chans')
    chans = request_details.chans;
    if ~iscell(chans)
        chans = {chans};
    end
else
    chans = {'BH?'};
end

% SNRmin = 5;

minmag = min(request_details.maglims);
maxmag = max(request_details.maglims);

gc_lims = request_details.gclims;     % [[4.7748;4.9821],[9.4758;9.8629]]  
              
redo_gc_lims = false;

ifsave = true;
if ~ifsave, datafile=''; end

chanstr = chans{1};
for ic = 2:length(chans)
    chanstr = [chanstr,',',chans{ic}];
end

javaaddpath('/Users/zeilon/Documents/MATLAB/IRIS-WS-2.0.15.jar')
addpath('~/Dropbox/MATLAB/seis_tools/breqfasting/')

%% get station info
stainfo = irisFetch.Stations('channel',network,station,'*',chanstr);
startdate = stainfo(1).StartDate;
enddate = stainfo(end).EndDate;

%% get event info
evinfo = irisFetch.Events('startTime',startdate,'endTime',enddate,...
        'radialcoordinates',[stainfo(1).Latitude,stainfo(1).Longitude,max(gc_lims), min(gc_lims)],...
        'maximumMagnitude',maxmag,'minimumMagnitude',minmag);
norids = length(evinfo);
evinfo_orig = evinfo;
    
%% cycle through events
% calc. rayp and expected arrival time for each 
gcarcs = zeros(norids,1);
seazs = zeros(norids,1);
arr_times = zeros(norids,length(phases));
rayps = zeros(norids,length(phases));
for ie = 1:norids
    [gcarcs(ie),seazs(ie)] = distance(stainfo(1).Latitude,stainfo(1).Longitude,...
                                evinfo(ie).PreferredLatitude,evinfo(ie).PreferredLongitude);
    for ip = 1:length(phases)
        tt = tauptime('z',evinfo(ie).PreferredDepth,'deg',gcarcs(ie),'ph',phases{ip});
        arr_times(ie,ip) = tt(1).time;
        rayps(ie,ip) = tt(1).rayparameter;
    end
end

%% limit by ray parm
if redo_gc_lims 
    figure(5),clf,set(gcf,'pos',[680   678   660   300])
    plot(gcarcs,seazs,'ob')
    title('P(b) S(r) - pick gc window','fontsize',20)
    yn = 'n';
    while ~strcmp(yn,'y')
    a = ginput(2); gc_lims = a(:,1); 
    yn = input(sprintf('%s gcarc limits: %.2f to %.2f. y/n  ',phases{ip},a(1,1),a(2,1)),'s');
    end
end

gdevts = find(gcarcs>=gc_lims(1) & gcarcs<=gc_lims(2));

%% good events only
evinfo = evinfo_orig(gdevts);
norids = length(evinfo);
arr_times = arr_times(gdevts,:);
gcarcs = gcarcs(gdevts);
seazs = seazs(gdevts);
rayps = rayps(gdevts,:);

% plot arrivals
figure(6), clf, hold on
for ip = 1:length(phases)
    if strcmp(phases{ip},'P')
        polar(d2r(seazs),rayp2inc(rayps(:,ip),8.1,6371-100),'ob');
    elseif strcmp(phases{ip},'S')
        polar(d2r(seazs),rayp2inc(rayps(:,ip),4.3,6371-100),'or');
    end
end


%% get data
SNR = zeros(norids,3,length(phases));
dat_all = zeros(diff(datwind)*samprate,norids,3,length(phases)); % channels in order ZRT
tt_all= zeros(diff(datwind)*samprate,norids,length(phases));
fprintf('Requesting %.0f events\n',length(evinfo))
for ie = 1:length(evinfo)
    fprintf('\nDownloading data for orid %.0f ',ie);
    for ip = 1:length(phases)
        fprintf('%s ',phases{ip})
        time0 = datestr(datenum(evinfo(ie).PreferredTime) + arr_times(ie,ip)/86400 + datwind(1)/86400 - 1/86400);
        time1 = datestr(datenum(evinfo(ie).PreferredTime) + arr_times(ie,ip)/86400 + datwind(2)/86400 + 1/86400);
    
%         tr_ev = irisFetch.Traces(network,station,'*','BH?',time0,time1,'includePZ');
        tr_ev = irisFetch.Traces(network,station,'*',chanstr,time0,time1); % grab PZ from resp file
        
        fprintf('...processing... ')
        
%        % continue if not 3 chans
%        if length(tr)<3, continue; end
%        if length(tr)>3 % two sensors/locations for this instrument
%            if mod(length(tr),3)~=0, continue; end % not multiple of 3 chans
%            while length(tr)>3
%                tr = tr([tr.depth] == min(unique([tr.depth]))); % prefer non-borehole
%                tr = tr([tr.sampleRate] == max(unique([tr.sampleRate]))); % prefer high samprate
%            end
%        end

        
        % continue if empty
        if isempty(tr_ev), fprintf('no data '); continue; end
        ischan = false(size(tr_ev));
        for ic = 1:length(tr_ev),ischan(ic)=~isempty(tr_ev(ic).network); end
        tr_ev = tr_ev(ischan); 
        if length(tr_ev)<3, fprintf('no data '); continue; end
        
        % find out how many different locations per station for this event
        [locs,nlocs] = tr_locs(tr_ev); 
        
        nlocav = 0; dat = [];
        for il = 1:nlocs % loop over locations
            tr = tr_ev(locs==il);
            
            % if 6 chans, and just repeat, kill three
            if length(tr)==6 
                if isequal(tr(1).data,tr(4).data), 
                    tr = tr(1:3);
                end
            end
            
            % continue if not 3 chans
            if length(tr)<3, fprintf('only %.0f chans ',length(tr)); continue; end
%             if mod(length(tr),3)~=0, continue; end % not multiple of 3 chans
            if length(tr)>3 % same location but still >3 traces - must be different samprates
                 tr = tr([tr.sampleRate] == max(unique([tr.sampleRate]))); % prefer high samprate
            end
            if length(tr)~=3, fprintf('%.0f chans ',length(tr)); continue; end
            

            %% down/up sample
			for ic = 1:3
				if samprate~=tr(1).sampleRate
					tr(ic).data = downsamp(tr(ic).data,tr(ic).sampleRate,samprate);
					tr(ic).sampleCount = length(tr(ic).data);
					tr(ic).sampleRate = samprate;
				end  
			end
		
			%% interp to common time axis - gets rid of requested buffer second
			tt_win = datwind(1) + [0:(diff(datwind)*tr(1).sampleRate - 1)]'./tr(1).sampleRate;
			datc = zeros(length(tt_win),3);
			for ic = 1:3
				tt_raw = tr(ic).startTime + [0:tr(ic).sampleCount-1]'./tr(ic).sampleRate/86400;
				tt_rel = (tt_raw - datenum(evinfo(ie).PreferredTime))*86400 - arr_times(ie,ip);
				datc(:,ic) = interp1(tt_rel,tr(ic).data,tt_win);
			end
		
			%% remove instrument response
			for ic = 1:3
                [pp,zz,gain] = resp_pz_getfromRESP(station,network,tr(ic).channel,tr(ic).location,tr(ic).startTime,respdir);
				datc(:,ic) = rm_resp(datc(:,ic),zz,pp,gain,tr(ic).sampleRate,0,0);
			end
		
			%% Rotate to ZNE 
			% rotation matrix for channels - Z(positive-DOWN), N, E directions 
			dp = [tr.dip];
			az = [tr.azimuth];
			Rmat = [sind(dp(1))  cosd(dp(1))*cosd(az(1))  cosd(dp(1))*sind(az(1));
					sind(dp(2))  cosd(dp(2))*cosd(az(2))  cosd(dp(2))*sind(az(2));
					sind(dp(3))  cosd(dp(3))*cosd(az(3))  cosd(dp(3))*sind(az(3))];

			datZNE = datc*Rmat;
			
			%% ZNE to ZRT
			dat = zne2zrt(datZNE,seazs(ie));% in ZRT, with T positive to right and Z positive down (still)
        
%         %% Processing
%         % detrend
%         dat = detrend(dat);
%         
%         % highpass filter at 100 s
%         dat = filt_quick(dat,1/100,samprate/2,1./samprate,2);
% 
%         %window
%         dat = flat_hanning_win(tt_win,dat,datwind(1),datwind(2),tapertime);

% %         compare
%         figure(45);clf;for ic = 1:3, hold on;plot(tr(ic).data);end
%         figure(44);clf; plot(tt_win,dat)
       
            nlocav = nlocav + 1;
        end % loop on locs
        
        if isempty(dat), continue; end
        %% average data across locations
        dat = sum(dat,3)./nlocav; 

        %% put in data structure
        dat_all(:,ie,:,ip) = dat;
        tt_all(:,ie,ip) = tt_win;
        fprintf(' stored')
    end % loop on events
end % loop on phases


%% kill nans or empties
gdevts = zeros(norids,length(phases));
for ie = 1:length(evinfo)
for ip = 1:length(phases)
    % check data not nan, and not all zeros
    gdevts(ie,ip) = all(all(~isnan(dat_all(:,ie,:,ip)))) & all(~all(dat_all(:,ie,:,ip)==0));
end
end
% good events only!
gdevts = find(all(gdevts,2));
evinfo = evinfo(gdevts);
norids = length(evinfo);
arr_times = arr_times(gdevts,:);
gcarcs = gcarcs(gdevts);
seazs = seazs(gdevts);
rayps = rayps(gdevts,:);
dat_all = dat_all(:,gdevts,:,:);
tt_all = tt_all(:,gdevts,:);

if isempty(tt_all), fprintf('>>>> No data\n');return; end

figure(66);plot(tt_all(:,1,1),sum(dat_all(:,:,1,1),2))

% %% align with CC
% for ip = 1:length(phases)
%     switch phases{ip}
%         case 'P', chid = 1; % align on Z
%         case 'S', chid = 3; % align on T (avoid R noise)
%     end
%     gdars = SNR(:,chid,ip)>SNRmin;
%     xcordat = filt_quick(dat_all(:,gdars,chid,ip),xcor_filtfs(1,ip),xcor_filtfs(2,ip),1./samprate);
%     xcortt = tt_all(:,1,ip);
% %     figure;plot(xcordat)
%     for ie = 1:sum(gdars),
%         dd = xcordat(:,ie);
%         ddmax = max(abs(dd));
%         ddsgn = sign( dd( find(abs(dd)==ddmax,1,'first') ) ); 
%         xcordat(:,ie) = dd.*ddsgn; % flip based on sign
%     end
%     figure(11)
%     plot(xcortt,xcordat);
%     xcordat = flat_hanning_win(tt_all(:,1,ip),xcordat,xcor_win(1),xcor_win(2),1);
%     figure(12)
%     plot(xcortt,xcordat);
%     [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat, 1./samprate, -datwind(1), 5,1)
%     
% % STACK raw traces to get overall arrtime
%     stk = zeros(length(xcortt),1);
%     for ig = 1:sum(gdars)
%         round(dcor(ig)*samprate);
%         % can only do this trick because padded with zeros
%         datstk = interp1(xcortt,xcordat(:,ig),xcortt+dcor(ig));
%         stk = stk + datstk(:);
%     end
%     figure(13), plot(xcortt,stk)
% %     pause
% end
arr_datenum_abs = zeros(norids,length(phases));
for ie = 1:norids
	arr_datenum_abs(ie,:) = arr_times(ie,:)/86400 + datenum(evinfo(ie).PreferredTime);
end

% figure(11)
% plot(tt_all(:,1,1),sum(dat_all(:,gdars,1,1),2))

%% save
eqar = struct('sta',station,'nwk',network,'phases',{phases},'components',{{'Z','R','T'}},...
              'slat',stainfo(1).Latitude,'slon',stainfo(1).Longitude,'selev',stainfo(1).Elevation,...
              'norids',norids,'elats',[evinfo.PreferredLatitude]','elons',[evinfo.PreferredLongitude]',...
              'edeps',[evinfo.PreferredDepth]','emags',[evinfo.PreferredMagnitude]','evtimes',{{evinfo.PreferredTime}'},...
              'gcarcs',gcarcs,'seazs',seazs,'rayps',rayps,'arr_times',arr_times,'arr_datenum_abs',arr_datenum_abs,...
              'dataZRT',dat_all,'tt',tt_all);
          
if ifsave
    fprintf('\nSAVING\n')
    datafile = sprintf('dat_%s_%s_%.0fto%.0f',station,network,gc_lims(1),gc_lims(2));
    save(datafile,'eqar');
end


    


end % on function
