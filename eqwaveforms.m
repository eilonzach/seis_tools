dbdir='/Volumes/Portadb/';
dbnam='portapng';
starttime=1; % time from orid to begin window
endtime=1600; % time from orid to end window
chan='BHZ';
filtsecs=[10 50];

orid=222;
sta='AGAN';

%% db details
db=dbopen(strcat(dbdir,dbnam),'r');
dbwf=dblookup_table(db,'wfdisc');
dbor=dblookup_table(db,'origin');
dbsite=dblookup_table(db,'site');
%% orid details
% orids=dbgetv(dbor,'orid')
% orid=input('Choose orid: ');
dbor.record=dbfind(dbor,sprintf('orid==%u',orid));
edep  =dbgetv(dbor,'depth');
elat  =dbgetv(dbor,'lat');  
elon  =dbgetv(dbor,'lon');
evtime=dbgetv(dbor,'time');
t0=evtime+starttime;
t1=evtime+endtime;
%% sta details
% stas=dbgetv(dbsite,'sta')
% sta=input('Choose sta: ','s');
dbsite.record=dbfind(dbsite,sprintf('sta=="%s"',sta));
slat  =dbgetv(dbsite,'lat');
slon  =dbgetv(dbsite,'lon');
%% arrivals
[delta, AZ] = distance(slat,slon,elat,elon);
[times, phasenames] = arrtimes( delta, edep );
fprintf('Phase %s arrives at %.0f\n',phasenames{1},times(1));
%% get data
[tt, dat, chans, nsamps, samprate, wfids] = dbgetwfz(db, sta, t0, t1);

filtfreqs=sort(1./filtsecs);
dt=unique(mean(diff(tt)));
tt=tt-evtime;

%% process waveforms
% FILTER %%% filt filt?
[b,a]=butter(3,filtfreqs.*2.*dt);
data_e=filtfilt(b,a,dat(:,1));
data_n=filtfilt(b,a,dat(:,2));
data_z=filtfilt(b,a,dat(:,3));
amp2=max([max(abs(data_e)),max(abs(data_n)),max(abs(data_z))]);
data_e=data_e./amp2;
data_n=data_n./amp2;
data_z=data_z./amp2;
tt0=tt;

% WINDOW - taper
we=window(@tukeywin,length(data_e));
wn=window(@tukeywin,length(data_n));
wz=window(@tukeywin,length(data_z));
data_e=data_e.*we;
data_n=data_n.*wn;
data_z=data_z.*wz;

%% plot data
figure(5)
plot(tt,data_z,'k','LineWidth',1.4)
[xlims]=ginput(2);
xlim(xlims(:,1));
ylim([-1.1 1.1]);
