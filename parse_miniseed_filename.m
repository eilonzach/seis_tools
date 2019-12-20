function [mseedfile_info] = parse_miniseed_filename(mseedfilename)
%[miniseedfile_info] = parse_miniseed_filename(mseedfilename)
%   
% Function to parse a miniseed file name into its constituent parts

if ~iscell(mseedfilename)
    mseedfilename = {mseedfilename}; 
end
for ifile = 1:length(mseedfilename)
        
%% Assume mseedfilename name format is '%NW.STAT.centaur-#_CTSN_YYYYMMDD_HHMMSS.miniseed'
C = textscan(mseedfilename{ifile},'%2s.%4s.centaur-%1u_%4u_%8s_%6s.miniseed');
try
nwk = C{1}{1};
sta = C{2}{1};
centaur_nchan = C{3};
centaur_SN = C{4};
yr = str2num(C{5}{1}(1:4));
mm = str2num(C{5}{1}(5:6));
dd = str2num(C{5}{1}(7:8));
HH = str2num(C{6}{1}(1:2));
MM = str2num(C{6}{1}(3:4));
SS = str2num(C{6}{1}(5:6));
jday = doy(yr,mm,dd);
% [nwk,strrem] = strtok(mseedfilename{ifile},'.');
% [sta,strrem] = strtok(strrem,'.');
% [chan,strrem] = strtok(strrem,'.');
% [loc,strrem] = strtok(strrem,'.');
% [yrstr,strrem] = strtok(strrem,'.');
% [jdaystr,strrem] = strtok(strrem,'.');
% [hhmmssstr] = strtok(strrem,'.');



starttime = datenum([yr,mm,dd,HH,MM,SS]);
mseedfile_info(ifile,1) = struct('sta',sta,'nwk',nwk,...
                      'year',yr,'jday',jday,'month',mm,'day',dd,...
                      'hour',HH,'min',MM,'sec',SS,'serialstart',starttime,...
                      'centaur_nchan',centaur_nchan,'centaur_SN',centaur_SN);
catch
    fprintf('! Problem parsing %s\n',mseedfilename{ifile});
end
end

end

