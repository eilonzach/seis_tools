% Script to remove data from folders if the waveforms associated with the
% files have been excised from the wfdisc in the antelope database - scans
% through data, checks if it is in the wfdisc and if not removes the files
% to a null directory

dbnam='portapng';
dbdir='/Volumes/zeilon/PNG/portapng/';
datadir='/Volumes/zeilon/PNG/portapng/data/';
nulldir='/Volumes/zeilon/PNG/portapng/data/dead/';
bufftime=60; % seconds after orid to check if file exists

allyall=0;

D=dir(datadir);
for id=4:length(D) % loop over stations
%% sort out sta names
if strcmp(D(id).name,'.')==1 ||...
            strcmp(D(id).name,'..')==1 ||...
                    strcmp(D(id).name,'.DS_Store')==1 % discount '.' files
continue
end
if ~isempty(regexp(D(id).name,'rot', 'once'))==1;
    sta=strtok(D(id).name,'rot');
    fprintf('STATION: %s (rotated)\n',sta);
elseif ~isempty(regexp(D(id).name,'rot', 'once'))==0;
    sta=D(id).name;
    fprintf('STATION: %s\n',sta);
end

%% get file names and details
%override all
if allyall==0,      yn='nall';
elseif allyall==1,  yn='yall';
end

F=dir(strcat(datadir,D(id).name));
for ifl=3:length(F) % loop over events, get all event & chan details
file=F(ifl).name; fprintf('\nFile: %s... ',file');
fullfile=strcat(datadir,D(id).name,F(ifl).name);

[out]=textscan(file,strcat('%4u.%3u.%2u.%2u.%2u.',sprintf('%s',sta),'.sac.%s'));
yyyy=out{1}; jjj=out{2}; hh=out{3}; mm=out{4}; ss=out{5}; chanx=char(out{6});
evtime=str2epoch(sprintf('%4u/%3u %2u:%2u:%2u',yyyy,jjj,hh,mm,ss));
if      strcmp(chanx,'e')==1, chan='BHE';
elseif  strcmp(chanx,'n')==1, chan='BHN';
elseif  strcmp(chanx,'z')==1, chan='BHZ';
end

wfsubstr1=sprintf('sta=="%s" && chan=="%s"',sta,chan);
wfsubstr2=sprintf('time < %.0f && endtime > %.0f',[1 1]*(evtime+bufftime));

db=dbopen(strcat(dbdir,dbnam),'r');
dbwf=dblookup_table(db,'wfdisc');
dbj1=dbsubset(dbwf,wfsubstr1);
dbj2=dbsubset(dbj1,wfsubstr2);
nrec=dbquery(dbj2,'dbRECORD_COUNT'); fprintf('nrec = %u... ',nrec);
if nrec>0
wfid=dbgetv(dbj2,'wfid');
end
dbclose(db)

% delete old file and wfdisc row query?
if nrec==0 
    if strcmp(yn,'yall')==0
	yn=input('Do you want to delete file? (y/n/yall) ','s');
    end
    if strcmp(yn,'y')==1 || strcmp(yn,'yall')==1 
    
	fprintf('excising file')
	movefile(fullfile,nulldir);
    
    elseif strcmp(yn,'n')==1
        continue
    end
end

end %loop on files
end %loop on stas


