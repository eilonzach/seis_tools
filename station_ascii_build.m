detail = 'station'; % 'station', 'channel', or 'response
nw = '2A';
stas = '*';
location = '*';
channel = 'BH*';

starttime = '1996 06 01'; % 'YYYY-MM-DD'
endtime = '2003 06 01'; % 'YYYY-MM-DD'
labo = -90;
lato = 90;
lole = -180;
lori = 180;

odir = '~/Work/Proposals/PacificArray/plotting/'; % need final slash
datfile = 'GLIMPSE_Pacificstas.txt';


%% ============ NO NEED TO CHANGE ANYTHING BELOW HERE =============== 
stadat  = irisFetch.Stations(detail,nw,stas,location,channel);

% parse for conditions
gd = true(length(stadat),1);
for is = 1:length(stadat)
    if str2epoch(stadat(is).EndDate)<str2epoch(starttime),
        gd(is) = false; continue
    end
    if str2epoch(stadat(is).StartDate)>str2epoch(endtime),
        gd(is) = false; continue
    end    
    if ~(stadat(is).Latitude >= labo && stadat(is).Latitude <= lato)
        gd(is) = false; continue
    end
    if lole<lori
        if ~(stadat(is).Longitude >= lole && stadat(is).Longitude <= lori)
            gd(is) = false; continue
        end
    elseif lole>lori
        if (stadat(is).Longitude>lole && stadat(is).Longitude<lori) || (stadat(is).Longitude<lole && stadat(is).Longitude>lori) 
            gd(is) = false; continue
        end
    end
end

stadat = stadat(gd)';

% figure; plot([stadat.Longitude],[stadat.Latitude],'o')

%% print to out file
fid = fopen([odir,datfile],'w+');
fprintf(fid,'sta   nwk   slat   slon   selev   startdate    enddate    status\n');
for is = 1:length(stadat)
    fprintf(fid,'%-6s,%2s,%8.4f,%9.4f,%7.1f,%10s,%10s,%6s\n',...
            stadat(is).StationCode,stadat(is).NetworkCode,...
            stadat(is).Latitude,stadat(is).Longitude,stadat(is).Elevation,...
            epoch2str(str2epoch(stadat(is).StartDate),'%Y-%m-%d'),...
            epoch2str(str2epoch(stadat(is).EndDate),'%Y-%m-%d'),...
            stadat(is).RestrictedStatus);
end
