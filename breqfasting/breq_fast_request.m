function breq_fast_request(label,usrname,stas,chans,nwks,locs,starttimes,endtimes,datatype,ofile)
%     breq_fast_request(label,stas,chans,nwks,locs,starttimes,endtimes,datatype,ofile)
% function to make a BREQ_FAST request file and email it to IRIS to request
% the data, by specifying a station and event parameters. This will create
% a request file with each data window requested on a separate line, of the
% format:
%   stat nw yyyy mm dd HH MM SS.sss  yyyy mm dd HH MM SS.sss #chan chans loc
% 
%
% INPUTS:
%  label   - identifying label for request (string)
%            IF empty, default is USER_request_TIMESTRING
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
%  endtimes - end time of requested data window. This can be given as a
%           [Nw x 1] vector of serial times (days since 0000-00-00) for Nw
%           different time windows. It can also be given as an {Nw x 1}
%           cell array of timestrings 
%           (suggested format 'yyyy-mm-dd HH:MM:SS.ssss')
%  datatype - type of data to request. Options:
%           'SEED', 'dataless_SEED', 'miniSEED', 'sync'
%            IF empty, default is 'SEED'
%  ofile -  name of output file containing the text of the request. 
%            IF empty, default is '~/Documents/MATLAB/label_BREQFAST_REQUEST'
% 
% If Ns==Nw, the script will assume each time goes with a single station,
% and the request will be Nrec=Ns(=Nw) lines long. Otherwise, the script
% will loop over times and stations to make a request with Nrec=NsxNw lines


%% fill in gaps
if isempty(nwks)
    nwks = {'*'};
end

if isempty(usrname)
    usrname = 'Researcher';
end

if isempty(label)
    [~,usr] = system('echo $USER'); usr = usr(~isspace(usr));
    label = [usr,'_request_',datestr(now,'yyyymmddHHMM')];
end
if isempty(locs)
    locs = '';
end

if nargin < 9 || isempty(datatype)
    % data type to request - if no match, sends to me
    % NB IF YOU DON'T WANT TO SEND AN EMAIL, WRITE DATATYPE = 'NULL'
    datatype = 'SEED'; % options are {'SEED','dataless SEED','miniSEED','sync'}
end

if nargin < 10 || isempty(ofile)
    ofile=['~/Documents/MATLAB/',label,'_BREQFAST_REQUEST'];
end

addpath([fileparts(which('breq_fast_process')),'/matguts']);

%% make sure in right formats
if ischar(stas), stas = {stas}; end
if ischar(nwks), nwks = {nwks}; end
if ~ischar(locs), locs = sprintf('%02s',locs); end
if ischar(locs), locs = {locs}; end
if ischar(starttimes), starttimes = {starttimes}; end
if ischar(endtimes), endtimes = {endtimes}; end
% convert time-string datetimes to serial time
starttimes = datenum(starttimes);
endtimes = datenum(endtimes);

%% request details
usrname(isspace(usrname)) = '_';
institution = 'UC Santa Barbara';
snail_mail  = '2116 Webb Hall, UCSB, Santa Barbara, CA 93106';
e_mail      = getpref('Internet','E_mail');
workphone   = '(805) 893-4688';
workfax     = '(805) 893-2314';

%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% SHOULDN'T NEED TO CHANGE ANYTHING BELOW HERE
%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Nsta = length(stas);
Ntws = length(starttimes);
Nnws = length(nwks);
Nlcs = length(locs);
if Nsta>1 && Nnws>1 && Nsta~=Nnws
    error('No. of stations and time windows not compatible')
end

%% header information in file
fid=fopen(ofile,'w+');
fprintf(fid,'.NAME %s\n',usrname);
fprintf(fid,'.INST %s\n',institution);
fprintf(fid,'.MAIL %s\n',snail_mail);
fprintf(fid,'.EMAIL %s\n',e_mail);
fprintf(fid,'.PHONE %s\n',workphone);
fprintf(fid,'.FAX %s\n',workfax);
fprintf(fid,'.MEDIA: Electronic\n');
fprintf(fid,'.ALTERNATE MEDIA: EXABYTE\n');
fprintf(fid,'.ALTERNATE MEDIA: DVD\n');
fprintf(fid,'.LABEL %s\n',label);
fprintf(fid,'.QUALITY B\n');
fprintf(fid,'.END\n\n');

%% For each event, run through stations getting time windows
if Nsta~=Ntws % loop over both stations and time windows
    Nreq = Nsta*Ntws;
    for ista=1:Nsta
    for itws=1:Ntws
        if Nnws==1; inws = 1; else, inws=ista; end
        if Nlcs==1; ilcs = 1; else, ilcs=ista; end
        % write to file
        fprintf(fid,'%s ',stas{ista}); % sta
        fprintf(fid,'%s ',nwks{inws});
        fprintf(fid,'%s  ',datestr(starttimes(itws),'yyyy mm dd HH MM SS.FFF'));
        fprintf(fid,'%s ',datestr(endtimes(itws),'yyyy mm dd HH MM SS.FFF'));
        fprintf(fid,'%u ',length({chans}));
        fprintf(fid,'%s ',chans);
        fprintf(fid,'%s\n',locs{ilcs});
    end
    end
elseif Nsta == Ntws % assume each station goes with corresponding time window
    Nreq = Nsta;
    for ista=1:Nsta
        if Nnws==1; inws = 1; else, inws=ista; end
        if Nlcs==1; ilcs = 1; else, ilcs=ista; end
        % write to file
        fprintf(fid,'%s ',stas{ista}); % sta
        fprintf(fid,'%s ',nwks{inws});
        fprintf(fid,'%s  ',datestr(starttimes(ista),'yyyy mm dd HH MM SS.FFF'));
        fprintf(fid,'%s ',datestr(endtimes(ista),'yyyy mm dd HH MM SS.FFF'));
        fprintf(fid,'%u ',length({chans}));
        fprintf(fid,'%s ',chans);
        fprintf(fid,'%s\n',locs{ilcs});
    end
end




%% read back into cell structure
nlines = 13 + Nreq;
messagetext=cell(nlines,1);

frewind(fid)
for il = 1:nlines
    messagetext(il,1) = eval('{fgets(fid)}');
end
   
fclose(fid);
%% 'to' email address
if strcmp(datatype,'SEED')==1; toaddress='breq_fast@iris.washington.edu'; 
elseif strcmp(datatype,'dataless_SEED')==1; toaddress='dataless@iris.washington.edu'; 
elseif strcmp(datatype,'miniSEED')==1; toaddress='miniseed@iris.washington.edu'; 
elseif strcmp(datatype,'sync')==1; toaddress='sync@iris.washington.edu'; 
else toaddress=e_mail; 
end
%% send to IRIS
sendmail(toaddress,'IRIS data request',messagetext)
sendmail(e_mail,'IRIS data request',messagetext)