function breq_fast_dataless_PZprocess( label,usrname,respdir,chans,ifdelete )
% breq_fast_dataless_PZprocess( label,usrname,respdir,ifdelete )
% function to download the output of a BREQ_FAST request that has
% been accepted and has a dataless SEED files waiting on an ftp server.  
% This script downloads that dataless seed and extracts the PZ response
% files.
%
% INPUTS:
%  label   - identifying label for request (string)
%            IF empty will fail
%  usrname - string identifying the user (e.g. 'Zach Eilon'). Spaces will
%            be replaced by underscores. 
%  respdir - directory in which to copy all the Response files
%  chans   - cell array of channel codes to keep
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
if isempty(usrname)
    usrname = 'Researcher';
end
usrname(isspace(usrname)) = '_';
if nargin < 3 || isempty(respdir)
    respdir = [wd,'/RESP'];
    mkdir(respdir);
end
if nargin < 4 || isempty(chans)
    chans = {'*'};
end
if nargin < 5 || isempty(ifdelete)
    ifdelete = true;
end

%% Get data from server
fprintf('   connecting to server\n')
ftpobj = ftp('ftp.iris.washington.edu','anonymous','eilon@ucsb.edu');
cd(ftpobj,['pub/userdata/',usrname]);
sf=struct(ftpobj);  sf.jobject.enterLocalPassiveMode();
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
    fprintf('SEED file not on server yet\n'); 
    traces = [];
    return
end

fprintf('   downloading %s\n',SEED_file_name)
file_name = mget(ftpobj,SEED_file_name);
close(ftpobj)

%% read the data with rdseed
for ii = 1:length(file_name)
    system(sprintf('/usr/local/bin/rdseed -f %s -p',file_name{ii}))
end

%% move to respdir
respfiles = dir;
cpresp = false(length(respfiles),1);
for ii = 1:length(respfiles)
	if ~any(regexp(respfiles(ii).name,'SAC_PZs')), continue; end
    fprintf('  parsing %s\n',respfiles(ii).name)
    for ic = 1:length(chans)
        if ~any(regexp(respfiles(ii).name,regexptranslate('wildcard',chans{ic}))), continue; end
        cpresp(ii) = true;
    end
end
respfiles = respfiles(cpresp);
for ii = 1:length(respfiles)
    fprintf('  copying %s to respdir\n',respfiles(ii).name)
    system(sprintf('cp %s %s',respfiles(ii).name,respdir));
end


%% clean up
cd(wd);
pause(0.1);
if ifdelete
    rmdir('temp_dir','s');
end

end

