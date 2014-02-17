% Builds a -.wfmeas file for the antelope "trexcerpt" function, so as to
% excerpt waveforms in explicit mode according to the rubric:
%     Use  the  times from table trxin.txt
%     to excerpt explicit waveform segments from  the  dbwf  database.
%     Each  row  in trxin.txt specifies a station, channel, and start time
%     for the excerpt.
%Inputs: a HarvardCMT.mat file, with cmt structure
%        a station file with rows e.g. (copy precisely; tabs/spaces matter)
%       turnontime      drift       latitude    longitude   depth   name
%     yyyy:jjj:hh:mm:ss	#.######	##.#### 	###.####	####	NAME 
% writes to outfile

%%CORRECTION will be used with awk for the simpler trexcerpt convert mode 
%output will be relecant station, channel, starttime, 

infile='/Users/Zach/Documents/MATLAB/portapng/statinfo.txt';
outfile='/Users/Zach/Documents/MATLAB/portapng/trxin_portapng3.txt';
% phases={'SKS,SKKS'};
runup = 100; %number of seconds before predicted SKS arrival to start window
rundown = 100; %number of seconds after predicted SKKS arrival to finish window
Mwmin = 6.0; %minimum magnitude
ondate = [ 2010 03 01]; %datevector of beginning of time window [ yyyy mm dd ] 
offdate = [ 2011 07 01 ]; %datevector of end of time window [ yyyy mm dd ] 
minang = 91;
maxang = 180;
%time after event to window
windstart = 0;
windend = 1800;



%%  Get event information from HarvardCMT.mat file
load('harvardCMT.mat'); %loads the cmt structure
goodev=find(cmt.Mw > Mwmin); %impose magnitude constraint
infoevs={'year','month','day','jjj','hour','minute','sec','lat','long','depth'};
for j=1:length(infoevs)
    assignin('base',char(infoevs(j)),eval(sprintf('cmt.%s(goodev)',char(infoevs(j)))));
    eval(sprintf('ev%s=%s;',char(infoevs(j)),char(infoevs(j))));
end

%% Get station information
fid=fopen(infile,'r');
% [time,drift,lat,long,depth,#name]...
D=textscan(fid,'%s %f %f %f %f %s','delimiter','\t');
fclose(fid);
stations=D{6}';
infostats={'lat','long','depth'}; infocolumns=[3,4,5];
for j=1:length(infostats)
    assignin('base',sprintf('sta%s',char(infostats(j))),D{infocolumns(j)});
end

chans={'BH0','BHE','BH1','BHN','BHZ'};
% chans={'BH0'}; % temporarily

% NB need to get land station information - ultimately need to get out of
% the .txt environment and start using the antelope-matlab
% interface

fout = fopen(outfile,'w'); % or 'a' for append
% Header
fprintf(fout, '%s \t %s\t %s\t %s  \t %s\n','sta','chan','yr-jday  hr:mi:sec','windl','phase');
fclose(fout);

nop=length(phases);
nos=length(stations);
noe=length(evdepth);
noc=length(chans);

for ip=1:nop
for is=1:nos
for ie=1:noe

    [foraz,backaz,gcdist,gcarc] = ellipsep(evlat(ie),evlong(ie),stalat(is),stalong(is));
    if  gcarc > minang...
     && gcarc < maxang...
     && datenum(evyear(ie),evmonth(ie),evday(ie)) > datenum(ondate)...
%% To window by certain phase arrival...
%      && datenum(evyear(ie),evmonth(ie),evday(ie)) < datenum(offdate)
%                 
%     tt=taupTime([],evdepth(ie),phases{ip},...
%         'sta',[stalat(is) stalong(is)],...
%         'evt',[evlat(ie) evlong(ie)]);
% 
%     skstime=tt(1).time;
%     skkstime=tt(2).time;
%     windstart=skstime-runup;
%     windend=skkstime+rundown;
%% No nonsense approach. Should always capture all useful phases
    winddur=windend-windstart;
    
%% parse starttime into yyy-jjj hh:mm:ss
    ws=rem(windstart,60);
    wmi=floor(windstart/60);
    
    SEC=rem(evsec(ie)+ws,60); MIpl=floor((evsec(ie)+ws)/60);
    MI=rem(evminute(ie)+wmi+MIpl,60); HRpl=floor((evminute(ie)+wmi+MIpl)/60);
    HR=rem(evhour(ie)+HRpl,24); DAYpl=floor((evhour(ie)+HRpl)/24);
    DAY=evday(ie)+DAYpl;
%     Jday=(datenum(eval(num2str(evyear(ie))),eval(num2str(evmonth(ie))),eval(num2str(DAY)))-datenum(eval(num2str(evyear(ie))),0,0));
    Jday=evjjj(ie)+DAYpl;
    
    % parse window length into hh:mm:ss
    wls=ceil(rem(winddur,60));
    wlmi=floor(rem(winddur,3600)/60);
    wlh=floor(winddur/3600);
    

    %% annoying things to get all time strings into right formats: 
    %annoying things to get starttime d.o.y. in right format (3-digit string)
    if Jday >= 100
        assignin('base','Jday',sprintf('%s',num2str(Jday)));
    elseif Jday >= 10 && Jday < 100
        assignin('base','Jday',sprintf('0%s',num2str(Jday)));
    elseif Jday < 10
        assignin('base','Jday',sprintf('00%s',num2str(Jday)));
    end
    %annoying things to get starttime HR in right format (2-digit string)
    if HR >= 10
        assignin('base','HR',sprintf('%s',num2str(HR)));
    elseif HR < 10
        assignin('base','HR',sprintf('0%s',num2str(HR)));
    end
    %annoying things to get starttime MI in right format (2-digit string)
    if MI >= 10
        assignin('base','MI',sprintf('%s',num2str(MI)));
    elseif MI < 10
        assignin('base','MI',sprintf('0%s',num2str(MI)));
    end
    %annoying things to get starttime SEC in right format (2-digit string)
    SEC=round(SEC);
    if SEC >= 10
        assignin('base','SEC',sprintf('%s',num2str(SEC)));
    elseif SEC < 10
        assignin('base','SEC',sprintf('0%s',num2str(SEC)));
    end
    %annoying things to get windl wlh in right format (2-digit string)
    if wlh >= 10
        assignin('base','wlh',sprintf('%s',num2str(wlh)));
    elseif wlh < 10
        assignin('base','wlh',sprintf('0%s',num2str(wlh)));
    end
    %annoying things to get windl wlmi in right format (2-digit string)
    if wlmi >= 10
        assignin('base','wlmi',sprintf('%s',num2str(wlmi)));
    elseif wlmi < 10
        assignin('base','wlmi',sprintf('0%s',num2str(wlmi)));
    end
    %annoying things to get windl wls in right format (2-digit string)
    if wls >= 10
        assignin('base','wls',sprintf('%s',num2str(wls)));
    elseif wls < 10
        assignin('base','wls',sprintf('0%s',num2str(wls)));
    end

    
 %% write to outfile
for ic=1:noc
fout = fopen(outfile,'a'); % or 'a' for append
% line of output
fprintf(fout,'%s\t %s\t %s\t %s\t %s\n',char(stations(is)),char(chans(ic)),sprintf('%s-%s %s:%s:%s',num2str(evyear(ie)),Jday,HR,MI,SEC),...
    sprintf('%s:%s:%s',wlh,wlmi,wls),...
    char(phases(ip)));
fclose(fout);
end;
end;
end;
end;
end;





%% Bear in mind the trexconv4splt build script
% #!/bin/bash
% # Updated 11/1/2012 10:12
% # This script requires a file (trxin.txt) of stations, channels, and start times for the waveform window 
% ## statA chanA timeA
% ## statB chanB timeB
% ## statC chanC timeC
% ## etc. 
% 
% # stat* is the name of the station (must match the antelope tables)
% # chan* is the name of the channel (must match the antelope tables)
% # time* is the time for the beginning of the window around the relevant arrival,
% #       in the format YYYY-DDD HH:MM:SS
% # The trexcerpt command builds sac files from the window accordingly
% 
% ## Required inputs - the name of the database, the channel names,
% ## and the file with the correction angles
% db=/Volumes/zeilon/PNG/cdpapuall	  #for example
% chanE=BH0	#for example
% chanN=BH1	#for example
% chanZ=BHZ	#for example
% echo "Give the name of the file containing the time excerpt data"
% read INFILE
% i=1
% ## go through line by line getting information
% while [ $i -le  $(awk 'END { print NR }' < $INFILE) ]; do
% stat=$(awk -v ii=$i 'NR==ii {print $1}'< $INFILE)
% chan=$(awk -v ii=$i 'NR==ii {print $2}'< $INFILE)
% timestart=$(awk -v ii=$i 'NR==ii {printf("%s %s",$3,$4)}'< $INFILE)
% windl=$(awk -v ii=$i 'NR==ii {print $5}'< $INFILE)
% 
% if [ "$chan" = "$chanE" ]; then
% CHAN=e
% elif [ "$chan" = "$chanN" ]; then
% CHAN=n
% elif [ "$chan" = "$chanZ" ]; then
% CHAN=z
% fi
% 
% # do the excerpting
% trexcerpt -vvo sc -w EXCERPT/%{sta}/%Y.%j.%H.%M.%S.%{sta}.sac.$CHAN\
%  -c "chan=='$chan' && sta=='$stat'" $db trialwf "$timestart" $windl
% # N.B. If a verbose output is not required, remove the "-v" in the two lines above
% # N.B. If confirmations for each change are not required, remove the "-c" in the two lines above
% let i=i+1
% done