function [ model ] = read_cardfile( cardfile )
% [ model ] = read_modelfile( cardfile )
fid = fopen(cardfile,'r');
modname = textscan(fid,'%s\n',1);
mdeets = fgetl(fid);
discs = fgetl(fid);
[nlay] = disc_vals(discs);
MOD = textscan(fid,'%8.0f %8.2f %8.2f %7.2f %7.1f %7.1f %8.2f %7.2f %.4f %*[^\n]',nlay);
fclose(fid);

model = struct( 'R',flipud(MOD{1}/1e3),...
                'rho',flipud(MOD{2}),...
                'Vp',flipud(MOD{3}),...
                'Vs',flipud(MOD{4}),...
                'Qk',flipud(MOD{5}),...
                'Qm',flipud(MOD{6})); 
model.eta = zeros(size(model.R));
model.Z = 6371 - model.R;

model.details = model_details(mdeets);
[model.nlay,model.discz] = disc_vals(discs,flipud(model.Z));




end


function mdeets = model_details(mdeets)
% function to parse model details from cardfile header
parsed_deets = sscanf(mdeets,'%f');
mdeets = struct('ifanis',parsed_deets(1)==1,...
                'cf',    parsed_deets(2),...
                'a',     parsed_deets(3)==1,...
                'b',     parsed_deets(3)==1);
end


function [nlay,zd] = disc_vals(discs,Z)
% function to parse discontinuity depths from cardfile header
zvals = sscanf(discs,'%u');
nlay = zvals(1);
    if nargin>1
        zvals = zvals(zvals<nlay);
        zd = Z(zvals);
    else
        zd=[];
    end
end