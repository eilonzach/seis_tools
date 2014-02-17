% Script to make a RECORDHEADERS format file for input into fortran
% programs like QnmA


% station specifications
chans={'BHE','BHN','BHZ'}; % e.g. BHE BHN BHZ
nchans=3;
samprate=50; % nominal samprate - will try to get this from antelope tables
network='ZN';
dostas={'DIOD'}; % default is {} - then will do all stations. Else, only does stations specified
ignorestas={'BBU1','KIR1','ESA1','SEHA','MGO','MISA','PMG'... 
            'SAGI','VAKU','WANI'}; % list of stations to ignore

% data type to request - if no match, sends to me
% NB IF YOU DON'T WANT TO SEND AN EMAIL, WRITE DATATYPE = 'NULL'
datatype = 'SEED'; % options are {'SEED','dataless SEED','miniSEED','sync'}

% Values to get time window
phase='S';
pretime=30; % seconds before predicted phase arrival to start window
posttime=60; % seconds after predicted phase arrival to end window

% outfile
odir='/Users/Zach/Documents/MATLAB/RECORDFILES/';
ofile_pref='directS'; % prefix for outfile - rest of name dictated by event

% database information
dbdir='/Volumes/zeilon/PNG/portapng/';
dbnam='portapng';

% event criteria
MSmin = 6.5;                % magnitude min.
Mbmin = 0;                  % magnitude min.
startdate = [2010,03,01];   % in form [yyyy,mm,dd]
enddate   = [2011,07,01]; 	% in form [yyyy,mm,dd]
mindepth = 200;             % in km
maxdepth = 1000;            % in km
refpt  = [-10,150];         % reference point for distance constraints
mindeg = 30;                % min. angular distance
maxdeg = 60;                % max. angular distance    
minaz = 0;                  % backaz range
maxaz = 360;                % backaz range


%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% SHOULDN'T NEED TO CHANGE ANYTHING BELOW HERE
%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% Get good events
db=dbopen(strcat(dbdir,dbnam),'r');
[goodorids]=gdorids_fn(db,MSmin,Mbmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz);
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
    % make new file
    fid=fopen(strcat(odir,ofile_pref,'.',epoch2str(evtime,'%Y_%j_%H_%M')),'w+');
    fprintf('Making file for orid %u, %s\n',orid,epoch2str(evtime,'%Y_%j_%H_%M'))
for is=1:nstas
    sta=char(stas(is));
    
    % if dostas specified, only do these
    if isempty(dostas)~=1
        if any(strcmp(sta,dostas))~=1
        continue
        end
    end
    % don't bother for unwanted stations
    if any(strcmp(sta,ignorestas))==1
        continue
    end

    
    % Window around arrival time
    tt=taupTime([],edep,phase,'sta',[slats(is),slons(is)],'evt',[elat,elon]);
    arrtime=tt.time;
    
    t0=evtime+arrtime-pretime;
    t1=evtime+arrtime+posttime;
    tt=t1-t0;
    
for ic=1:nchans
    chan=char(chans(ic));
    if strcmp(chan(end),'E')==1; hang=90; vang=0;
    elseif strcmp(chan(end),'N')==1; hang=0; vang=0;
    elseif strcmp(chan(end),'Z')==1; hang=0; vang=-90;
    end
    
	db=dbopen(strcat(dbdir,dbnam),'r');
    dbsensor=dblookup_table(db,'sensor');
    dbinst=dblookup_table(db,'instrument');
    dbsensor.record=dbfind(dbsensor,sprintf('sta=="%s" && chan=="%s"',sta,chan));
    inst_no=dbgetv(dbsensor,'inid');
    dbinst.record=dbfind(dbinst,sprintf('inid==%u',inst_no));
    if dbquery(dbinst,'dbRECORD_COUNT')==1
    samprate=dbgetv(dbinst,'samprate');
    end
    dbclose(db)
    nsamps=tt*samprate + 1;
    
    % stn-nw, chan, slat, slon, selev, burial depth, hang, vang, samprate, nsamps, year, julday, hr, min, sec, total time, "class"

    % write to file
    fprintf(fid,sprintf('%4s',sta)); % sta
    fprintf(fid,'-%s  ',network); % nw
    fprintf(fid,'%s      ',chan); % nw
    fprintf(fid,'%8.4f ',slats(is)); % slat
    fprintf(fid,'%9.4f ',slons(is)); % slon
    fprintf(fid,'%5.0f.    0.0  ',selevs(is)*1000); % selev, burial depth
	fprintf(fid,'%4.1f  %5.1f     ',hang, vang); % nw
	fprintf(fid,'%8.4f ',samprate); % nw
    fprintf(fid,'%7.0f ',nsamps); % nw
    fprintf(fid,'%s  ',epoch2str(t0,'%Y %j %H %M %S.%s'));
    fprintf(fid,'%9.2f ',tt); % nw
    fprintf(fid,'Data\n');
end
end
    fclose(fid);
%     return % only do for first event
end