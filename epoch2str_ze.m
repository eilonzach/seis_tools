function [ str ] = epoch2str_ze( epochaltime,strformat )
% [ str ] = epoch2str_ze( epochaltime,strformat )
% 
% my function to do the epoch to string conversion, for situations where
% the matlab-antelope toolbox isn't working. Relies on the inbuilt datestr
% function in matlab.
% 
% Z. Eilon, April 2017

serialtime = epoch2serial_ze(epochaltime);

% get strformat from epoch2str format into datestr format
strformat = regexprep(strformat,'%A','dddd');
strformat = regexprep(strformat,'%a','ddd');
strformat = regexprep(strformat,'%B','mmmm');
strformat = regexprep(strformat,'%b','mmm');
strformat = regexprep(strformat,'%D','mmddyy');
strformat = regexprep(strformat,'%d','dd');
strformat = regexprep(strformat,'%Y','yyyy');
strformat = regexprep(strformat,'%y','yy');
strformat = regexprep(strformat,'%m','mm');
strformat = regexprep(strformat,'%H','HH');
strformat = regexprep(strformat,'%M','MM');
strformat = regexprep(strformat,'%S','SS');
strformat = regexprep(strformat,'%s','FFF');
strformat = regexprep(strformat,'%p','PM');

str = datestr(serialtime,strformat);

end

