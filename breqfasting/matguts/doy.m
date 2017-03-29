function [JJJ] = doy(YYYY,MM,DD,mode)
%[JJJ] = doy(YYYY,MM,DD,mode) or doy([YYYY,MM,DD],[],[],mode)
% converts calendar date (day + month) into day of year
% if mode is 's', will output a 3-character string
% if mode is 'i', will output a number
% 
% Written by Zach Eilon, 2011, updated 2017

if nargin==1 || (isempty(MM) && isempty(DD))
    DD = YYYY(3);
    MM = YYYY(2);
    YYYY = YYYY(1);
end

if ischar(YYYY), YYYY = str2double(YYYY);end
if ischar(MM), MM = str2double(MM);end
if ischar(DD), DD = str2double(DD);end

if nargin < 4 || isempty(mode)
    mode = []; 
end

JJJ=datenum(YYYY,MM,DD)-datenum(YYYY,0,0);

if strcmp(mode,'s')
	JJJ=sprintf('%03s',num2str(JJJ));
end

