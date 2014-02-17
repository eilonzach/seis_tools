function [strike,dip,rake] = cmtdetail(evtime,elat,elon)
%[strike,dip,rake] = cmtdetail(time,lat,long)
% script to pull out the important event details from cmt catalogue. Will
% attempt to find event that matches the time to within the day, the lat
% and lon to within a tenth of a degree
%
%INPUTS
% evtime - in epochal format
% elat
% elon

load('/Users/Zach/Documents/MATLAB/SplitLab1.0.5/harvardCMT.mat')
yes_lat=find(abs(cmt.lat - elat)<0.1); %impose lat constraint
yes_lon=find(abs(cmt.long - elon)<0.1); %impose lon constraint
time=datevec(epoch2str(evtime,'%Y/%m/%d %H:%M:%S')); %serial date
yes_time=find(abs(datenum(cmt.year,cmt.month,cmt.day,cmt.hour,cmt.minute,cmt.sec)...
    -datenum(time))<2);
goodev=intersect(yes_lat,yes_lon); goodev=intersect(goodev,yes_time);

if length(goodev)>1
fprintf('Multiple events satisfy criteria:\n')
fprintf('Input event: %s (%6.2fN, %6.2fE)\n',epoch2str(evtime,'%Y/%m/%d %H:%M:%S'),elat,elon)
fprintf('Found events:\n')
for ie=1:length(goodev)
    kk=goodev(ie);
    fprintf('%u)  %4u/%u/%u %2u:%2u:%2u  (%6.2fN, %6.2fE)  Mb=%3.1f depth=%.1fkm\n',...
        ie,cmt.year(kk),cmt.month(kk),cmt.day(kk),cmt.hour(kk),cmt.minute(kk),cmt.sec(kk),...
        cmt.lat(kk),cmt.long(kk),cmt.Mb(kk),cmt.depth(kk))
end
ans1=input('Enter number corresponding to desired event: ');
goodev=goodev(ans1);
end
if isempty(goodev)==1, fprintf('No events seem to match\t');
else
strike=cmt.strike(goodev);
dip=cmt.dip(goodev);
rake=cmt.rake(goodev);
end


end

