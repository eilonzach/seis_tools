function [ model ] = read_seismodel_nc( ncfile )
%[ model ] = read_seismodel_nc( ncfile )
%  
%  Function to use netcdf toolbox tools to read a netcdf ".nc" file and
%  extract a seismic model in a matlab-friendly structure.
% 
% INPUTS:
%  ncfile  - name of the .nc model file
% 
% OUTPUTS:
%  model  - structure with fields (where applicable):
%   .Z
%   .lat
%   .lon
%   .Vsv
%   .Anis
% 
% Z.E. 08/2016

model = struct('Z',[],'lat',[],'lon',[],'Vsv',[],'Vsh',[],'Vpv',[],'Vph',[],'rho',[]);

%% get ncid
ncid = netcdf.open(ncfile);

%% get variables
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
vars = cell(nvars,1);
for iv = 1:nvars
    vars{iv} = netcdf.inqVar(ncid,iv-1);
end

% fns = fieldnames(model);

for iv = 1:nvars
    switch vars{iv}
        case {'depth','Z','z'}
            model.Z = double(netcdf.getVar(ncid,iv-1));
        case {'latitude','lat'}
            model.lat = double(netcdf.getVar(ncid,iv-1));
        case {'longitude','lon'}
            model.lon = double(netcdf.getVar(ncid,iv-1));
        case {'Vsv','vsv','vs','Vs'}
            model.Vsv = double(netcdf.getVar(ncid,iv-1));
        case {'Vsh','vsh'}
            model.Vsh = double(netcdf.getVar(ncid,iv-1));  
        case {'Vpv','vpv','vp','Vp'}
            model.Vpv = double(netcdf.getVar(ncid,iv-1));
        case {'Vph','vph'}
            model.Vph = double(netcdf.getVar(ncid,iv-1));
        case {'rho','density'}
            model.rho = double(netcdf.getVar(ncid,iv-1));
    end
end

model.nz = length(model.Z);
model.nlat = length(model.lat);
model.nlon = length(model.lon);
model.npts = model.nz*model.nlat*model.nlon;
          

            
netcdf.close(ncid)
end

