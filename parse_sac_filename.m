function [sacfile_info] = parse_sac_filename(sacfilename)
%[sacfile_info] = parse_sac_filename(sacfilename
%   
% Function to parse a sac file name into its constituent parts

if ~iscell(sacfilename)
    sacfilename = {sacfilename}; 
end
for ifile = 1:length(sacfilename)
        
%% Assume SACfile name format is '%NW.STAT..CHX.D.YYYY.JJJ.HHMMSS.SAC'    
% C = textscan(sacfilename{ifile},'%2s.%4s..%3s.D.%4u.%3u.%6s.SAC');
% nwk = C{1}{1};
% sta = C{2}{1};
% chan = C{3}{1};
% yr = C{4};
% jday = C{5};
% HH = str2num(C{6}{1}(1:2));
% MM = str2num(C{6}{1}(3:4));
% SS = str2num(C{6}{1}(5:6));
[nwk,strrem] = strtok(sacfilename{ifile},'.');
[sta,strrem] = strtok(strrem,'.');
[chan,strrem] = strtok(strrem,'.');
[loc,strrem] = strtok(strrem,'.');
[yrstr,strrem] = strtok(strrem,'.');
[jdaystr,strrem] = strtok(strrem,'.');
[hhmmssstr] = strtok(strrem,'.');

yr = str2num(yrstr);
jday = str2num(jdaystr);
[mm,dd] = calday(yr,jday);
HH = str2num(hhmmssstr(1:2));
MM = str2num(hhmmssstr(3:4));
SS = str2num(hhmmssstr(5:6));

starttime = datenum([yr,mm,dd,HH,MM,SS]);
sacfile_info(ifile,1) = struct('sta',sta,'nwk',nwk,'chan',chan,'loc',loc,...
                      'year',yr,'jday',jday,'month',mm,'day',dd,...
                      'hour',HH,'min',MM,'sec',SS,'serialstart',starttime);

end

end

