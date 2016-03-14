function [ reloc ] = readhypodd( reloc_file )
%[ reloc ] = readhypodd( reloc_file )
%   simple funtion to read a hypodd .reloc file, producing a structure that
%   contains the data about the relocated earthquakes

fid = fopen(reloc_file,'r');
A = textscan(fid,'%n %f %f %f %f %f %f %f %f %f %n %n %n %n %n %f %f %n %n %n %n %f %f %n');
fclose(fid);
reloc = struct('orid',A{1},'lat',A{2},'lon',A{3},'depth',A{4},'lat_error',A{9},'lon_error',A{8},'depth_error',A{10},...
               'mag',A{17},'ncatP',A{20},'ncatS',A{21},'rmscat',A{23},...
               'datenum',datenum(A{11},A{12},A{13},A{14},A{15},A{16}),...
               'datestr',{cellstr(datestr(datenum(A{11},A{12},A{13},A{14},A{15},A{16})))},...
               'epochtime',str2epoch(cellstr(datestr(datenum(A{11},A{12},A{13},A{14},A{15},A{16}),31))),...
               'NCCP',A{18},'NCCS',A{19},'rmsCC',A{22},'clusteri',A{24});


end
