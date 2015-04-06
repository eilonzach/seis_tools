function [minA,x,y] = mingrid(A)
% [minA,x,y] = mingrid(A)
%MINGRID This script finds the (x,y) location and value of the minimum of a
%MxN matrix, e.g. an error grid where:
% x is the column number of the minimum - the x cooridinate
% y is the row number of the minimum - the y coordinate
if nargin==1
    [~,x]=min(min(A));
    [minA,y]=(min(A(:,x)));
end
end

