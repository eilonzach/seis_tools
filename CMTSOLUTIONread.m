function [eq] = CMTSOLUTIONread(CMTSOLUTIONfile)
% [eq] = CMTSOLUTIONread(CMTSOLUTIONfile)
%  function to read in a CMTsolution file and spit out all the useful
%  information - BEWARE times. The origin time of the event is, properly,
%  equal to evtime-halfdur.
%
% INPUTS
%  CMTSOLUTIONfile - string with name and full path of file
% 
% OUTPUTS
%  eq      - structure with all earthquake information
%  eq.time - [yr, mon, day, hr, min, sec] NB THIS IS CENTROID TIME
%  eq.lat  - elat (degrees)
%  eq.lon  - elon (degrees)
%  eq.dep  - edep (km)
%  eq.M    - 3x3 full moment tensor (dyne*cm) (r, theta, phi) (+Z,S,E)
%  eq.halfdur - half duration (sec) 
%  eq.location - name of earthquake region
%  eq.name - GCMT name
% 
%       Z. Eilon	21 May 2013

fid = fopen(CMTSOLUTIONfile);
eq = textscan(fid,'%4s %f%f%f%f%f%f %f %f %f %f%f %s\n',1,'delimiter','\n');
location = char(eq{13});
eq = cell2mat(eq(2:12));
evname = textscan(fid,'event name: %s\n',1,'delimiter','\n');
tshift = textscan(fid,'time shift: %f\n',1,'delimiter','\n');
halfdur = textscan(fid,'half duration: %f\n',1,'delimiter','\n');
elat = textscan(fid,'latitude: %f\n',1,'delimiter','\n');
elon = textscan(fid,'longitude: %f\n',1,'delimiter','\n');
edep = textscan(fid,'depth: %f\n',1,'delimiter','\n');
MM = textscan(fid,'%s%f\n',6);
fclose(fid);

[evtime] = datevec(datenum([eq(1:5) eq(6)+cell2mat(tshift)]));

M = [MM{2}(1) MM{2}(4) MM{2}(5)
     MM{2}(4) MM{2}(2) MM{2}(6)
     MM{2}(5) MM{2}(6) MM{2}(3)];

eq = struct([]);
eq(1).time = evtime;
eq.lat = cell2mat(elat);
eq.lon = cell2mat(elon);
eq.dep = cell2mat(edep);
eq.M = M;
eq.halfdur = cell2mat(halfdur);
eq.location = location;
eq.name = char(evname{1}); eq.name = eq.name(isspace(eq.name)~=1);


end
