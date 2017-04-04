function [ serialtime ] = epoch2serial_ze( epochaltime )
%  [ serialtime ] = epoch2serial_ze( epochaltime )
% 
% my function to do the epoch to serial time conversion, for situations where
% the matlab-antelope toolbox isn't working. Serial time is given in days
% since Jan 01 0000; epochal time in seconds since Jan 01 1970.
% 
% Z. Eilon, April 2017

stJ12970 = datenum(1970,1,1,0,0,0);

serialtime = epochaltime/24/3600 + stJ12970;

end


