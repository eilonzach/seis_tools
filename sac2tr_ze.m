function traces = sac2tr_ze(sac_files,xml_files)
% traces = sac2tr_ze(sac_files,xml_files)
%
% Function to turn a SAC file (or list of SAC files) into a MATLAB traces
% structure, akin to what was returned from iris_Fetch.Traces.
% If there are instrument response files in xml format associated with the
% SAC files, then those can be listed in a second input. This list can
% either be the same length as the sac_files list, with the .xml file in
% each list element matching the corresponding sac file, or the list can
% simply include all the stations and this script will search for the
% matching file. 

arguments
    % Nx1 cell array of strings/char-vecs, or Nx1 array of strings, or single character vector. 
    % Default is everything in a directory.
    % Does accept wildcards
    sac_files = "*.sac";
    % cell array of strings/char-vecs, or array of strings, or single character vector. Default is everything in a directory.
    xml_files = "";
end

addpath('~/Dropbox/MATLAB/seis_tools') % for readsac

%% parse sac_files
sac_files_input = sac_files;

% get the sac files into an Nx1 array of strings
% if just one character array, convert to string
if ischar(sac_files), sac_files = string(sac_files); end
% if cell array, convert to array of strings
if iscell(sac_files)
    sac_files_new = "";
    for iff = 1:length(sac_files)
        sac_files_new(iff,1) = string(sac_files{iff});
    end
    sac_files = sac_files_new;
end


%% parse xml_files
xml_files_input = xml_files;

% get the xml files into an Nx1 array of strings
% if just one character array, convert to string
if ischar(xml_files), xml_files = string(xml_files); end
% if cell array, convert to array of strings
if iscell(xml_files)
    xml_files_new = "";
    for iff = 1:length(xml_files)
        xml_files_new(iff,1) = string(xml_files{iff});
    end
    xml_files = xml_files_new;
end

if length(xml_files)==length(sac_files) % looks like there is a match
    if_xml_match_sac = true;
else 
    if_xml_match_sac = false;
end


%% loop over sac files and load
traces = [];
for iscf = 1:numel(sac_files)

    % first process files with dir
    sac_file_i = dir(sac_files(iscf));

    % % Read all SACs, keep only BHZ/BHN/BHE
    recs = [];
    % used_paths = strings(0,1);
    for k = 1:numel(sac_file_i)
        sac_path = fullfile(sac_file_i(k).folder, sac_file_i(k).name);
        R = readsac(sac_path);
        if isempty(R), warning('Skipping unreadable SAC: %s', sac_path); continue; end
        chan = upper(strtrim(R.KCMPNM));
        if numel(chan)==3 && (strncmp(chan,'BH',2) || strncmp(chan,'HH',2)) && any(chan(3) == ['Z','N','E','1','2'])
            R.path = sac_path;
            recs = [recs; R]; %#ok<AGROW>
        end
    end
    if isempty(recs)
        warning('No BHZ/BHN/BHE/BH1/BH2 or HHZ/HHN/HHE/HH1/HH2 or  channels found. Nothing to load.');
    end

    % parse out networks
    unwks = string(unique({recs.KNETWK}));
    for inwk = 1:length(unwks)
        recsn = recs(string({recs.KNETWK}) == unwks(inwk));
    
    % parse out stations
    ustas = string(unique({recsn.KSTNM}));
    for ista = 1:length(ustas)
        recsns = recsn(string({recsn.KSTNM}) == ustas(ista));


    % if both BH? and HH?, then just HH?
    havechans = string({recsns.KCMPNM});
    okchans = true(length(havechans),1);
    if any(ismember(havechans,"HHE")),okchans(ismember(havechans,"BHE"))=0;end
    if any(ismember(havechans,"HHN")),okchans(ismember(havechans,"BHN"))=0;end
    if any(ismember(havechans,"HHZ")),okchans(ismember(havechans,"BHZ"))=0;end
    if any(ismember(havechans,"HH1")),okchans(ismember(havechans,"BH1"))=0;end
    if any(ismember(havechans,"HH2")),okchans(ismember(havechans,"BH2"))=0;end
    recsns = recsns(okchans);

    % For each channel, if multiple files, choose highest sample rate, then longest duration
    chans = {'Z','N','E','1','2'};
    for ic = 1:numel(chans)
        C = chans{ic};
        rr = recsns(contains(havechans,C));
        if isempty(rr), continue; end
        % what if multiple files - must be split up traces
        if numel(rr) > 1
            continue
            samprate = 1./[rr.DELTA];
            dur = [rr.endTime] - [rr.startTime];
            [~, pick] = max(samprate*1e6 + dur);
            rr = rr(pick);
        end

        % rename data field if needed
        rr.data = rr.DATA1;
        rr = rmfield(rr,'DATA1');

        samprate = round(1./rr.DELTA);

        t0 = sac_start_datenum(rr);
        tend = t0 + (numel(rr.data)-1)/samprate/86400;
    

        tr = struct();
        tr.network  = rr.KNETWK;
        tr.station  = rr.KSTNM;
        tr.location = rr.KHOLE;
        if strcmp(tr.location,' '),tr.location = '';end
        tr.channel  = rr.KCMPNM;
        tr.quality  = '';
        
        % get associated xml file
        xml_file_i = [];
        nwstastr = sprintf('%2s_%s',tr.network,tr.station);
        if if_xml_match_sac % 1:1 matching of sac and xml files
            xml_file_i = dir(xml_files(iscf));
        else % see if there is a match to one of the listed xml files, use that
            try
                xml_file_i = dir(xml_files(contains(xml_files,nwstastr)));
            end
        end
        % otherwise, hunt for file in the sac directory
        if isempty(xml_file_i)
            xml_file_i = dir(fullfile(sac_file_i(k).folder,...
                [nwstastr,'*.xml']));
        end

        try 
            ifxml = true;
            xml_path = fullfile(xml_file_i.folder, xml_file_i.name);
        catch
            xml_path =  fullfile(sac_file_i(k).folder,...
                [nwstastr,'*.xml']);
            ifxml = false;
        end

        try
            SX = readstruct(xml_path,'FileType','xml');
        catch
            warning('readstruct failed on %s, skip', xml_path); 
            ifxml = false;
        end
    
        if ifxml
            chmeta = read_stationxml_meta(SX, tr.network, tr.station, tr.location, tr.channel, t0, tend);
            if isempty(chmeta) || isempty(fieldnames(chmeta)), ifxml = false; end
        end 

        if ifxml
            tr.latitude  = chmeta.Latitude.Text;
            tr.longitude = chmeta.Longitude.Text;
            tr.elevation = chmeta.Elevation;
            tr.depth     = chmeta.Depth;
            tr.azimuth   = chmeta.Azimuth.Text;
            tr.dip       = chmeta.Dip.Text;
        else
            tr.latitude  = rr.STLA;
            tr.longitude = rr.STLO;
            tr.elevation = rr.STEL;
            tr.depth     = rr.STDP;
            tr.azimuth   = nan;
            tr.dip       = nan;
        end


        if ifxml
            [sens,sfreq,su,instr] = extract_sensitivity_rs(chmeta);
            tr.sensitivity = sens; 
            tr.sensitivityFrequency = sfreq; 
            tr.sensitivityUnits = su; 
            tr.instrument = instr;
        else
            tr.sensitivity = nan; 
            tr.sensitivityFrequency = nan; 
            tr.sensitivityUnits = nan; 
            tr.instrument = rr.KINST;
        end

        tr.data = rr.data(:);
        tr.sampleCount = numel(rr.data);
        tr.sampleRate = samprate;
        tr.startTime = datenum(datestr(t0,'yyyy-mm-dd HH:MM:SS.FFF'));
        tr.endTime   = datenum(datestr(tend,  'yyyy-mm-dd HH:MM:SS.FFF'));

        if ifxml
            tr.sacpz = extract_sacpz_rs(chmeta);
        else
            tr.sacpz = [];
        end

        try
            traces = [traces; tr]; %#ok<AGROW>
        catch
        keyboard
        end
    end % loop on channels to try
    end % loop on unique stations in nwk
    end % loop on unique networks 


end % loop on sac files



end % function



%% subfunctions
function v = get_sac_str(H, name, def)
    v = def; cands = {name, lower(name), upper(name)};
    for i=1:numel(cands)
        fn = cands{i};
        if isfield(H,fn)
            x = H.(fn);
            if ischar(x) || isstring(x), v = strtrim(char(x)); return; end
            if isnumeric(x), v = strtrim(char(string(x))); return; end
        end
    end
end

function v = get_sac_num(H, name, def)
    v = def; cands = {name, lower(name), upper(name)};
    for i=1:numel(cands)
        fn = cands{i};
        if isfield(H,fn)
            x = H.(fn);
            if isnumeric(x), v = double(x); return; end
            if (ischar(x) || isstring(x)) && ~isempty(x)
                t = str2double(char(x)); if ~isnan(t), v = t; return; end
            end
        end
    end
end

function t0 = sac_start_datenum(H)
    yr = get_sac_num(H,'NZYEAR',NaN);
    jd = get_sac_num(H,'NZJDAY',NaN);
    hh = get_sac_num(H,'NZHOUR',0);
    mm = get_sac_num(H,'NZMIN',0);
    ss = get_sac_num(H,'NZSEC',0);
    ms = get_sac_num(H,'NZMSEC',0);
    if isnan(yr) || isnan(jd), t0 = NaN; return; end
    t0 = datenum(datetime(yr,1,1) + days(jd-1) + hours(hh) + minutes(mm) + seconds(ss) + milliseconds(ms));
end

function C = read_stationxml_meta(SX, net, sta, loc, cha, t0, tend)
    C = struct();
    % Navigate FDSNStationXML -> Network -> Station -> Channel
    nets = getlist(SX,'Network');
    for i=1:numel(nets)
        N = nets{i};
        if ~streq(gettext(N,'codeAttribute'), net), continue; end
        STS = getlist(N,'Station');
        for j=1:numel(STS)
            S = STS{j};
            if ~streq(gettext(S,'codeAttribute'), sta), continue; end
            CHS = getlist(S,'Channel');
            CHS = CHS{1};
            for k=1:numel(CHS)
                Ch = CHS(k);
                if ~streq(gettext(Ch,'codeAttribute'), cha), continue; end
                if ~streq(gettext(Ch,'locationCodeAttribute'), loc), 
                    if strcmp(loc,' ')&&isempty(gettext(Ch,'locationCodeAttribute'))
                        % fine, they are both flavours of empty
                    elseif eval(loc) == eval(gettext(Ch,'locationCodeAttribute'))
                        % fine, they have same numerical value, just
                        % different char number.
                    else
                        continue; % no, they differ
                    end
                end
                on = -inf; off = inf;
                sd = gettext(Ch,'startDateAttribute'); if ~isempty(sd), on = datenum(parse_iso(sd)); end
                ed = gettext(Ch,'endDate');   if ~isempty(ed), off = datenum(parse_iso(ed)); end
                if ~(tend <= on || t0 >= off)
                    % Attach station-level coords if missing
                    if ~isfield(Ch,'Latitude'),  Ch.Latitude  = getnum_struct(S,'Latitude'); end
                    if ~isfield(Ch,'Longitude'), Ch.Longitude = getnum_struct(S,'Longitude'); end
                    if ~isfield(Ch,'Elevation'), Ch.Elevation = getnum_struct(S,'Elevation'); end
                    C = Ch; return
                end
            end
        end
    end
end

function [sens,freq,units,instr] = extract_sensitivity_rs(C)
    sens = NaN; freq = NaN; units = ''; instr = '';
    if isempty(C) || ~isstruct(C), return; end
    % Instrument description
    instr = gettext(C,'description');
    if isempty(instr)
        sensFld = findfield(C,'Sensor');
        if ~isempty(sensFld) && isfield(C.(sensFld),'Description')
            instr = gettext(C.(sensFld),'Description');
        end
    end
    % InstrumentSensitivity
    respFld = findfield(C,'Response');
    if ~isempty(respFld) && isfield(C.(respFld),'InstrumentSensitivity')
        IS = C.(respFld).InstrumentSensitivity;
        sens = str2double(charify(getfield_flex(IS,'Value')));
        freq = str2double(charify(getfield_flex(IS,'Frequency')));
        try
            units = charify(getfield_flex(IS,'InputUnits'));
        catch
            units = charify(getfield_flex(getfield_flex(IS,'InputUnits'),'Name'));
        end
        if ~isnan(sens), return; end
    end
    % Fallback Stage[1]
    if ~isempty(respFld) && isfield(C.(respFld),'Stage')
        ST = C.(respFld).Stage; if ~iscell(ST), ST = {ST}; end
        if ~isempty(ST) && isfield(ST{1},'StageSensitivity')
            SS = ST{1}.StageSensitivity;
            sens = str2double(charify(getfield_flex(SS,'Value')));
            freq = str2double(charify(getfield_flex(SS,'Frequency')));
            units = charify(getfield_flex(SS,'InputUnits'));
        end
    end
end

function sacpz = extract_sacpz_rs(C)
    
    % sacpz = '';
    respFld = findfield(C,'Response');
    rsp = C.(respFld);
    if isempty(rsp) || ~isfield(rsp,'Stage'), return; end
    InsSens = rsp.InstrumentSensitivity.Value;
    for istg = 1:numel(rsp)
        if isfield(rsp.Stage(istg),'PolesZeros') && ~isempty(rsp.Stage(istg).PolesZeros)
             pz = rsp.Stage(istg).PolesZeros; break
        end
    end
    % get total gain
    GainProd = 1;
    for istg = 1:numel(rsp)
        GainProd = GainProd*rsp.Stage(istg).StageGain.Value;
    end
    % extract poles and zeros out of stupid fields
    pps = extract_poles(pz.Pole);
    zzs = extract_poles(pz.Zero);


    sacpz = struct('units',upper(pz.InputUnits.Name),...
                   'constant',pz.NormalizationFactor*GainProd,...
                   'poles',pps,'zeros',zzs);
            

end

function val = iff(cond,a,b), if cond, val=a; else, val=b; end, end

function v = getnum(S, name, default)
    v = default;
    if isempty(S) || ~isstruct(S), return; end
    fld = findfield(S, name);
    if ~isempty(fld)
        x = S.(fld);
        if isstruct(x) && isfield(x,'Text'), x = x.Text; end
        if ischar(x) || (isstring(x) && isscalar(x)), x = char(x); end
        if isnumeric(x), v = double(x);
        else, t = str2double(x); if ~isnan(t), v = t; end
        end
    end
end

function fld = findfield(S, want)
    % Case-insensitive, ignores namespace prefixes, matches end of tag
    fld = '';
    fns = fieldnames(S);
    for i=1:numel(fns)
        n = fns{i};
        base = n;
        % strip namespaces like fdsn:Latitude
        c = strfind(n,':'); if ~isempty(c), base = n(c(end)+1:end); end
        if strcmpi(base, want), fld = n; return; end
    end
end

function tf = streq(a,b)
    tf = strcmp(charify(a), charify(b));
end

function s = charify(x)
    if isstring(x) || ischar(x), s = char(x);
    elseif isstruct(x) && isfield(x,'Text'), s = char(x.Text);
    else, s = char(string(x));
    end
end

function L = getlist(S, tag)
    L = {};
    fld = findfield(S, tag);
    if isempty(fld) || ~isfield(S,fld), return; end
    x = S.(fld);
    if isempty(x), return; end
    if iscell(x), L = x; % seem to require x to be a cell - not wanted.
    else, L = {x};
    end
end

function t = gettext(S, tag)
    t = '';
    fld = findfield(S, tag);
    if isempty(fld) || ~isfield(S,fld), return; end
    x = S.(fld);
    if isstruct(x) && isfield(x,'Text'), t = char(x.Text);
    elseif ischar(x) || isstring(x), t = char(x);
    else, t = char(string(x));
    end
end

function v = getnum_struct(S, tag)
    v = NaN;
    fld = findfield(S, tag);
    if isempty(fld) || ~isfield(S,fld), return; end
    x = S.(fld);
    if isstruct(x) && isfield(x,'Text'), x = x.Text; end
    if isnumeric(x), v = double(x);
    else, v = str2double(char(string(x)));
    end
end

function dt = parse_iso(s)
    % Handles 'YYYY-MM-DDTHH:MM:SS[.fff][Z]'
    s = char(s);
    s = regexprep(s,'Z$','');
    try
        dt = datetime(s,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS','TimeZone','UTC');
    catch
        dt = datetime(s,'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC');
    end
end



function v = getfield_flex(S, tag)
    v = '';
    fld = findfield(S, tag);
    if isempty(fld) || ~isfield(S,fld), return; end
    x = S.(fld);
    if isstruct(x) && isfield(x,'Text'), v = x.Text; else, v = x; end
end



function p_or_z = extract_poles(p_or_z_struct)
    N = numel(p_or_z_struct);
    p_or_z = nan(N,1);
    for ip = 1:N
        p_or_z(ip) = getfield_flex(p_or_z_struct(ip),'Real') ...
                      + 1i*getfield_flex(p_or_z_struct(ip),'Imaginary');
    end
end
