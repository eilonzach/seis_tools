function [ norids,orids,elats,elons,edeps,evtimes,mags ] = dbor_asciiread( dbdir,dbnam )
%[ norids,orids,elats,elons,edeps,evtimes,mags ] = dbor_asciiread( dbdir,dbnam )
%   
% Function to bypass the Antelope-MATLAB interface (for computers that
% don't have it, or for versions of MATLAB that don't support it). This
% simple function directly reads the ascii file for the origin table of the
% specified database, and parses it into the useful vectors.

orfile = [dbdir,'/',dbnam,'.origin'];

fid = fopen(orfile,'r');
A = textscan(fid,'%8.4f %9.4f %9.4f %16.5f %8d %8d %8d %4d %4d %4d %8d %8d -  - %12.4f - %7.2f %8d %7.2f %8d %7.2f %8d -               -  %22d %16.5f\n');
fclose(fid);

orids = A{5};
norids = length(orids);
elats = A{1};
elons = A{2};
edeps = A{3};
evtimes = A{4};
mags = A{16};

end

