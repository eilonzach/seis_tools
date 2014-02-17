% Script to make a request file to send to breq_fast to get IRIS data, by
% specifying a station and event parameters - uses irisFetch tools
% ? If 'dataype' is set properly, will automatically send file to IRIS
% ? Uses all stations in antelope db (unless explicitly told to ignore)
% ? Uses all events that conform to set criteria
% requires site and origin tables

label='EilonPMGdata'; % also name of outfile

% station specifications
chans='BH?'; % e.g. BHE BHN BHZ or ? as a wildcard for any spot
nchans=3;
network='IU';
dostas={'ANMO'};
ignorestas={'BBU1','KIR1','ESA1','SEHA','MGO','MISA'}; % list of stations to ignore

% data type to request - if no match, sends to me
% NB IF YOU DON'T WANT TO SEND AN EMAIL, WRITE DATATYPE = 'NULL'
datatype = 'SEED'; % options are {'SEED','dataless SEED','miniSEED','sync'}

% Values to get time window
phase=''; % name of phase or LEAVE BLANK ('') to do evtime
pretime=120; % seconds before predicted phase arrival to start window
posttime=2400; % seconds after predicted phase arrival to end window

% outfile
odir='/Users/Zach/Documents/MATLAB/PNG_swsplit/IRIS PMG data/';
ofile=strcat(label,' IRIS REQUEST');

% database information
dbdir='/Users/Zach/Documents/MATLAB/';
dbnam=''; % or LEAVE BLANK ('') to do from refpt and CMT

% event criteria
MSmin = 0;                % magnitude min.
Mwmin = 6.5;                  % magnitude min.
startdate = [2010,04,01];   % in form [yyyy,mm,dd]
enddate   = [2012,07,31]; 	% in form [yyyy,mm,dd]
mindepth = 0;             % in km
maxdepth = 1000;            % in km
refpt  = [-10, 150];         % reference point for distance constraints
mindeg = 90;                % min. angular distance
maxdeg = 135;                % max. angular distance    
minaz = 0;                  % backaz range
maxaz = 360;                % backaz range

% request details
name        = 'Zachary Eilon';
institution = 'LDEO Columbia University';
snail_mail  = '61 Route 9W, Palisades, NY 10964';
e_mail      = 'zeilon@ldeo.columbia.edu';
workphone   = '(845) 365-8460';
workfax     = '(845) 365-8101';

%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% SHOULDN'T NEED TO CHANGE ANYTHING BELOW HERE
%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% Get good events

db=dbopen(strcat(dbdir,dbnam),'r');
[goodorids]=gdorids_fn(db,MSmin,Mwmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz);
dbclose(db)

norids=length(goodorids);

%% Get stations
db=dbopen(strcat(dbdir,dbnam),'r');
dbsite=dblookup_table(db,'site');
stas=dbgetv(dbsite,'sta');
slats=dbgetv(dbsite,'lat');
slons=dbgetv(dbsite,'lon');
selevs=dbgetv(dbsite,'elev');
nstas=dbquery(dbsite,'dbRECORD_COUNT');
dbclose(db)


%% header information in file
fid=fopen(strcat(odir,ofile),'w+');
fprintf(fid,'.NAME %s\n',name);
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
for ie=1:norids
    orid=goodorids(ie);
    db=dbopen(strcat(dbdir,dbnam),'r');
    dbor=dblookup_table(db,'origin');
    dbjor=dbsubset(dbor,sprintf('orid==%u',orid));
    
    evtime=dbgetv(dbjor,'time');
    elat=dbgetv(dbjor,'lat');
    elon=dbgetv(dbjor,'lon');
    edep=dbgetv(dbjor,'depth');
    dbclose(db)
    for is=1:nstas
        sta=char(stas(is));

        % don't bother for unwanted stations
        if any(strcmp(sta,ignorestas))==1
            continue
        end
        if isempty(dostas)~=1
            if any(strcmp(sta,dostas))~=1
                 continue
            end
        end

        if isempty(phase)~=1
        % Window around arrival time
        tt=taupTime([],edep,phase,'sta',[slats(is),slons(is)],'evt',[elat,elon]);
        arrtime=tt.time;
        else
            arrtime = 0;
        end

        t0=evtime+arrtime-pretime;
        t1=evtime+arrtime+posttime;
        % write to file
        fprintf(fid,sprintf('%s ',sta)); % sta
        fprintf(fid,'%s ',network);
        fprintf(fid,'%s  ',epoch2str(t0,'%Y %m %d %H %M %S.%s0'));
        fprintf(fid,'%s ',epoch2str(t1,'%Y %m %d %H %M %S.%s0'));
        fprintf(fid,'%u ',nchans);
        fprintf(fid,'%s\n',chans);
    end
end
%% read back into cell structure
nlines = 13 + nstas*norids;
messagetext=cell(nlines,1);

frewind(fid)
for il = 1:nlines
    messagetext(il,1) = eval('{fgets(fid)}');
end
   
fclose(fid);
%% 'to' email address
if strcmp(datatype,'SEED')==1; toaddress='breq_fast@iris.washington.edu'; 
elseif strcmp(datatype,'dataless SEED')==1; toaddress='dataless@iris.washington.edu'; 
elseif strcmp(datatype,'miniSEED')==1; toaddress='miniseed@iris.washington.edu'; 
elseif strcmp(datatype,'sync')==1; toaddress='sync@iris.washington.edu'; 
else toaddress=e_mail; 
end
%% send to IRIS
sendmail(toaddress,'IRIS data request',messagetext)
sendmail(e_mail,'IRIS data request',messagetext)