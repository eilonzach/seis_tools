function [iunq,unqstas,unqnwks] = stanw_unique(stas,nwks)
% [iunq,unqstas,unqnws] = stanw_unique(stas,nws)
% 
%   Given two [nsta,1] cell arrays of stations and networks that may have
%   repeats, parse out the unique combinations of both, and return an array
%   of indices of the unique combos

stanw = strcat(stas,nwks);
stanw = regexprep(stanw,' ','');
[~,iunq] = unique(stanw,'stable');

unqstas = stas(iunq);
unqnwks = nwks(iunq);


end

