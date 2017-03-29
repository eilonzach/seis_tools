function [ traces ] = breq_fast_process( label,usrname,stas,chans,nwks,locs,starttimes,ifdelete )
%[ traces ] = breq_fast_process( label,stas,chans,nwks,locs,starttimes )
% function to parse and download the output of a BREQ_FAST request that has
% been accepted and has SEED files waiting on an ftp server. This will
% download the data and output it as a [Nevt x Nchan] traces structure.
%
% INPUTS:
%  label   - identifying label for request (string)
%            IF empty will fail
%  usrname - string identifying the user (e.g. 'Zach Eilon'). Spaces will
%            be replaced by underscores. 
%  stas    - string or {Nsx1} cell array of station names
%  nwks    - string or {Nsx1} cell array of networks (must either be one
%            value or a cell array of the same length as the stas cell 
%            array (N.B. may be wildcard, '*')
%            IF empty, default is '*'
%  locs    - number, string or cell array of locations (must either be one
%            value or a cell array of the same length as the stas cell 
%            array (N.B. may be wildcard, '*')
%            IF empty, default is '00'
%  starttimes - beginning time of requested data. This can be given as a
%           [Nw x 1] vector of serial times (days since 0000-00-00) for Nw
%           different time windows. It can also be given as an {Nw x 1}
%           cell array of timestrings 
%           (suggested format 'yyyy-mm-dd HH:MM:SS.ssss')
%  ifdelete -  option to delete the SEED, SAC, and SAC_PZs files  
%            IF empty, default is true
% 
% If Ns==Nw, the script will assume each time goes with a single station,
% and the request will be Nrec=Ns(=Nw) lines long. Otherwise, the script
% will loop over times and stations to make a request with Nrec=NsxNw lines


wd = pwd;
addpath(wd);
mkdir('temp_dir')
cd('temp_dir')

%% fill in the blanks
if isempty(locs), locs = '*'; end
if isempty(nwks), nwks = '*'; end
if isempty(chans), chans = '*'; end
if isempty(usrname)
    usrname = 'Researcher';
end
if nargin < 8 || isempty(ifdelete)
    ifdelete = true;
end

usrname(isspace(usrname)) = '_';

%% Get data from server
fprintf('   connecting to server\n')
ftpobj = ftp('ftp.iris.washington.edu','anonymous','eilon@geol.ucsb.edu');
cd(ftpobj,['pub/userdata/',usrname]);
% parse to find SEED file that matches "label" and is most recent
aa = dir(ftpobj);
t0 = 0;
for ii=1:length(aa)
    if ~isempty(regexp(aa(ii).name,label,'once')) && aa(ii).datenum>t0
        SEED_file_name = aa(ii).name;
        t0 = aa(ii).datenum; % default to latest
    end
end

if ~exist('SEED_file_name','var')
    close(ftpobj)
    cd(wd);
    error('SEED file not on server yet'); 
end

fprintf('   downloading %s\n',SEED_file_name)
file_name = mget(ftpobj,SEED_file_name);
close(ftpobj)

%% read the data with rdseed
for ii = 1:length(file_name)
    system(sprintf('/usr/local/bin/rdseed -f %s -p',file_name{ii}))
    system(sprintf('/usr/local/bin/rdseed -f %s -d',file_name{ii}))
end


%% make sure in right formats
if ischar(stas), stas = {stas}; end
if ischar(nwks), nwks = {nwks}; end
if ~ischar(locs), locs = char(locs); end
if ischar(locs), locs = {locs}; end
if ischar(starttimes), starttimes = {starttimes}; end
if ischar(chans), chans = {chans}; end

Nsta = length(stas);
Ntws = length(starttimes);
Nnws = length(nwks);
Nlcs = length(locs);
if Nsta>1 && Nnws>1 && Nsta~=Nnws
    error('No. of stations and time windows not compatible')
end


if Nsta~=Ntws % loop over both stations and time windows
    Nreq = Nsta*Ntws;
elseif Nsta == Ntws % assume each station goes with corresponding time window
    Nreq = Nsta;
end

%% For each event, run through stations getting time windows
for ii=1:Nreq
    % complex sorting out of indices if vectors are different lengths;
    
    if Nsta == 1
        ista = 1; 
    else
        if Nsta == Ntws
            ista=ii;
        elseif Nsta~=Ntws
            ista = rem(ii-1,Nsta)+1;
        end   
    end
    
    if Ntws==1
        itws = 1; 
    else
        if Nsta == Ntws
            itws=ista;
        elseif Nsta~=Ntws
            itws = floor((ii-1)/Nsta)+1;
        end
    end
    
    if Nnws==1; inws = 1; else, inws=ista; end
    if Nlcs==1; ilcs = 1; else, ilcs=ista; end

    
    for ic = 1:length(chans)
        jday = doy(datestr(starttimes(itws),'yyyy'),datestr(starttimes(itws),'mm'),datestr(starttimes(itws),'dd'),'s');
        sactime = [datestr(starttimes(itws),'yyyy.'),jday,datestr(starttimes(itws),'.HH.MM.SS.*')];
        sacfile = [sactime,'.',nwks{inws},'.',stas{ista},'.',locs{ilcs},'.',regexprep(chans{ic},'?','*'),'.*.SAC'];
        chanfiles = dir(sacfile);
        for jj = 1:length(chanfiles)
            if exist('tr','var')==2, delete(tr); end
            tr = irisFetch.SAC2Trace(chanfiles(jj).name);
            tr.dip = tr.dip-90; % no idea why dips are 90? off...
            % now get the response information for this trace
            pzfile = sprintf('SAC_PZs_%s_%s_%s_%s*',tr.network,tr.station,tr.channel,tr.location);
            respfiles = dir(pzfile);
            for iresp = 1:length(respfiles)
                [~,~,~,~,~,rt0,rt1]=parse_sacpz_filename(respfiles(iresp).name);
                rt0 = datenum([rt0(1),jday2mody(rt0(1), rt0(2)),rt0(3),rt0(4),rt0(5)]);
                rt1 = datenum([rt1(1),jday2mody(rt1(1), rt1(2)),rt1(3),rt1(4),rt1(5)]);
                if (tr.startTime > rt0) && (tr.startTime < rt1), continue, end
            end
            if isempty(respfiles(iresp).name)
                stop
            end
            [zzs, pps, const, unit] = read_sac_pole_zero(respfiles(iresp).name);
            tr.sacpz.constant = const;
            tr.sacpz.poles = pps;
            tr.sacpz.zeros = zzs;
            tr.sacpz.units = unit;
            
            traces(ii,jj) = tr;
        end
    end
end
% add a null final trace structure if the final trace isn't there. 
if size(traces,1)~=Nreq
    nulltr = tr;
    fns = fieldnames(tr);
    for ifn = 1:length(fns), nulltr.(fns{ifn}) = []; end
    traces(Nreq,:) = nulltr;
end

cd(wd);
if ifdelete
rmdir('temp_dir','s');

end

end

