function radpplot_all(str,dip,rak,fign)
% [ Z, X, Y ] = radpplot(str,dip,rak[,fign]) 
% This script plots all three radiation patterns for the input fault 
% geometry
% 
%INPUTS
% str   = fault strike, degrees from north
% dip   = fault dip, degrees from horizontal
% rak   = fault rake, +90 is pure thrust, 0 is left lateral for a ss
% fign  = figure number to do plots in 
% subp  = subplot indices - vector [y x i]
%
% for plotting more on same axes, use conversion:
% r = r_max * (sin(inc/2)/sin(45))
% x = sin(az)*r';
% y = cos(az)*r';
% for this script, as is, rmax is 4
%
%   Z. Eilon    20 May 2013

if nargin<4
    fign=1;
end

for ii= 1:3
    radpplot(str,dip,rak,ii,fign,[1,3,ii]);
end
    set(gcf,'Position',[100 550 1500 400]);
end