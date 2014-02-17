%% Change formats of filenames
% Script to change the formats of the filenames in a given folder from the
% standard output when you run mseed2sac on IRIS data in mseed format
% This script changes it to the format splitlab likes:
% yyyy.jjj.hh.mm.ss.sac.x
fdir = '~/Documents/MATLAB/PNG_swsplit/IRIS PMG data/';
ichans = {'BHE','BHN','BHZ','BH1','BH2'};
ochans = 'enzne';

fprintf('CAREFUL WITH BH1/2/E/N !!!\N');
D = dir(fdir);
for ff = 3:length(D);
    ifile = D(ff).name;
    if length(ifile)<30 % skip if not right file format
        continue
    else
        %loop over parts of filename
        [~,rmn] = strtok(ifile,'.'); 
        [sta,rmn] = strtok(rmn,'.'); 
        if strcmp(rmn(2),'.')~=1; [~,rmn] = strtok(rmn,'.'); end % do extra if needed here
        [chan,rmn] = strtok(rmn,'.');   
            ochan = ochans(strcmp(chan,ichans));
        [~,rmn] = strtok(rmn,'.');
        [yr,rmn] = strtok(rmn,'.,'); 
        [jday,rmn] = strtok(rmn,','); 
            for ii = 1:3-length(jday); jday = strcat('0',jday); end
        [hr,rmn] = strtok(rmn,',:'); 
            for ii = 1:2-length(hr); hr = strcat('0',hr); end
        [min,rmn] = strtok(rmn,':'); 
            for ii = 1:2-length(min); min = strcat('0',min); end
        [sec,rmn] = strtok(rmn,':.'); 
            for ii = 1:2-length(sec); sec = strcat('0',sec); end
        if strcmp(rmn,'.SAC')~=1
            error('filename seems to be wrong format for file %u\n',ff)
        end
    end
    ofile = sprintf('%4s.%3s.%2s.%2s.%2s.%s.sac.%s',yr,jday,hr,min,sec,sta,ochan);
    movefile(strcat(fdir,ifile),strcat(fdir,ofile));
end
