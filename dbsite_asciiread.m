function [ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype,network ] = dbsite_asciiread( dbdir,dbnam )
%[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype,network ] = dbsite_asciiread( dbdir,dbnam )
%   
% Function to bypass the Antelope-MATLAB interface (for computers that
% don't have it, or for versions of MATLAB that don't support it). This
% simple function directly reads the ascii file for the site table of the
% specified database, and parses it into the useful vectors.
% 
% Z. Eilon, April 2017

sitefile = [dbdir,'/',dbnam,'.site'];

% count lines
nlines = 0;
fid = fopen(sitefile,'r');
while 1
    out = fgetl(fid); 
    if out==-1, break; end
    nlines = nlines+1;
end
fclose(fid);

nstas = nlines;
stas = cell(nstas,1);
slats = zeros(nstas,1);
slons = zeros(nstas,1);
selevs = zeros(nstas,1);
ondate = zeros(nstas,1);
offdate = zeros(nstas,1);
staname = cell(nstas,1);
statype = cell(nstas,1);
network = cell(nstas,1);

fid = fopen(sitefile,'r');
for is = 1:nstas
    txtl = fgetl(fid);
    
    stas{is} = sscanf(txtl(1:6),'%6s',1);
    
    B = sscanf(txtl(7:24),'%8d %8d',2);
    ondate(is) = B(1);
    offdate(is) = B(2);
    
    C = sscanf(txtl(26:54),'%9f %9f %9f',3);
    slats(is) = C(1);
    slons(is) = C(2);
    selevs(is) = C(3);
    
    staname{is} = txtl(56:105);
    statype{is} = sscanf(txtl(106:110),'%s');
    network{is} = sscanf(txtl(112:113),'%s');
end
fclose(fid);


end
