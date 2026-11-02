function [sacfile_info] = parse_sac_filename(sacfilename)
%[sacfile_info] = parse_sac_filename(sacfilename)
%   
% Function to parse a sac file name into its constituent parts
% sacfilename can either be a character array, or a cell with multiple
% character arrays.

if ~iscell(sacfilename)
    sacfilename = {sacfilename}; 
end
for ifile = 1:length(sacfilename)
        
%% Assume SACfile name format is '%NW.STAT.LOC.CHX.D.YYYY.JJJ.HHMMSS.SAC'    
% % FIRST WAY 
% C = textscan(sacfilename{ifile},'%2s.%4s..%3s.D.%4u.%3u.%6s.SAC');
% nwk = C{1}{1};
% sta = C{2}{1};
% chan = C{3}{1};
% yr = C{4};
% jday = C{5};
% HH = str2num(C{6}{1}(1:2));
% MM = str2num(C{6}{1}(3:4));
% SS = str2num(C{6}{1}(5:6));

% % SECOND WAY 
% [nwk,strrem] = strtok(sacfilename{ifile},'.');
% [sta,strrem] = strtok(strrem,'.');
% [loc,strrem] = strtok(strrem,'.');
% [chan,strrem] = strtok(strrem,'.');
% [qual,strrem] = strtok(strrem,'.');
% [yrstr,strrem] = strtok(strrem,'.');
% [jdaystr,strrem] = strtok(strrem,'.');
% [hhmmssstr] = strtok(strrem,'.');

% % THIRD WAY 
ttt = textscan(sacfilename{ifile},'%s','Delimiter','.'); ttt = ttt{1};
nwk = ttt{1};
sta = ttt{2};
loc = ttt{3};
chan = ttt{4};
qual = ttt{5};
yrstr = ttt{6};
jdaystr = ttt{7};
hhmmssstr = ttt{8};


yr = str2num(yrstr);
jday = str2num(jdaystr);
[mm,dd] = calday(yr,jday);
HH = str2num(hhmmssstr(1:2));
MM = str2num(hhmmssstr(3:4));
SS = str2num(hhmmssstr(5:6));

starttime = datenum([yr,mm,dd,HH,MM,SS]);
sacfile_info(ifile,1) = struct('sta',sta,'nwk',nwk,'chan',chan,'loc',loc,'quality',qual,...
                      'year',yr,'jday',jday,'month',mm,'day',dd,...
                      'hour',HH,'min',MM,'sec',SS,'serialstart',starttime);

end

end

