email_matlab_setup
addpath('matguts');
figure(1);clf; set(gcf,'pos',[97 369 1090 500]);
%%



%% set up parameters for the data request
label = 'breq_test';
name = 'Seis Mologist'; % no punctuation, spaces are allowed;
stas = {'KIEV','HRV'};
nwks = 'IU';
chans = 'BHZ';
locs = '00';
datatype = 'SEED';
ofile = 'request_details';
odir = '.';

% get data for Tohoku-Oki 2011 earthquake
starttimes = {'11 Mar 2011 05:00:00'};
endtimes = {'11 Mar 2011 07:00:00'};

%% requesting the data
breq_fast_request(label,name,stas,chans,nwks,locs,starttimes,endtimes,datatype,ofile)

%% pausing to give the DMC processing time
fprintf('NOW IRIS 5 minutes to respond... check your email!\n')
pause(5*60);
fprintf('Once IRIS emails to say the data is ready, run the\nrest of the script.\n')

%% Downloading and processing the data
[ tr ] = breq_fast_process( label,name,stas,chans,nwks,locs,starttimes,true );

%% giving data time vectors
for itr = 1:length(tr)
    if isempty(tr(itr).network), continue; end
    tr(itr).tt = [0:tr(itr).sampleCount-1]'/tr(itr).sampleRate/86400 + tr(itr).startTime;
end

%% plotting
figure(1)
plot(tr(1).tt,tr(1).data);
datetick('x',13)
set(gca,'fontsize',14,'xlim',[tr(1).tt(5e4),tr(1).tt(end)])
xlabel('Time (11 March 2011, UTC)','fontsize',17,'fontweight','bold');
ylabel('Counts','fontsize',17,'fontweight','bold')
title(sprintf('%s, %s',tr(1).station,tr(1).channel),'fontsize',18)
