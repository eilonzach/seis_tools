function [pp,zz,gain,ounit] = resp_pz_getfromRESP(sta,nwk,chan,loc,time_serial,respdir)
% [pp,zz,gain,ounit] = resp_pz_getfromRESP(sta,nwk,chan,loc,time_serial,respdir)
% 
%   Function to search a RESP directory for known station, event, channel
%   and event timeand get the appropriate SACPZ response file, read it in,
%   and spit out the vectors of poles, zeros, and the gain.

respfile_pre = ['SAC_PZs_',nwk,'_',sta,'_',chan,'_',loc];
respfile = dir([respdir,respfile_pre,'*']);
respfile = flipud(respfile); % use more recent starts first
if length(respfile) > 1 % multiple response files for this nwk+sta+chan
    for ir = 1:length(respfile)
    % get beginning and end of response file in serial time
    [~,~,~,~,~,resp_begin,resp_end]=parse_sacpz_filename(respfile(ir).name);
    [rbm,rbd] = calday(resp_begin(1),resp_begin(2));
    [rem,red] = calday(resp_end(1),resp_end(2));
    rb = datenum([resp_begin(1),rbm,rbd,resp_begin(3:end)]);
    re = datenum([resp_end(1),rem,red,resp_end(3:end)]);
    if (rb < time_serial) && (re > time_serial)
        % if event within this time, this is desired respfile
        break;
    end
    end
elseif length(respfile)==1
    ir = 1;
else
    fprintf('NO RESP'), return
end

% use function from seizmo toolbox
[pz]=readsacpz_rdseed([respdir,respfile(ir).name]);
gain = pz.k;
zz = complex(pz.z{1});
pp = complex(pz.p{1});
ounit = pz.input{1};
if strcmp(ounit,'M') % PZ file set up to convert to displacement
    % correct to VELOCITY by removing a zero
    zz = zz(2:end);
end            


end

