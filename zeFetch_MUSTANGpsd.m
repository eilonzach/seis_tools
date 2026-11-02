function [psdstruct] = zeFetch_MUSTANGpsd(nwk,sta,loc,chan,time_start,time_end,ifplot)
% [psdstruct] = zeFetch_MUSTANGpsd(nwk,sta,loc,chan,time_start,time_end,ifplot)
% 
% Function to get all individual hourly spectra for a station for some date
% range. Note, will always request quality 'M', unless you go into code to
% change.
%
% INPUTS
%  nwk  = two-character code
%  sta  = four(or five) character code
%  loc = '*' % two char code
%  chan = '?H?'; % three char code
%  time_start = '2010-02-27T06:30:00.000' % string
%  time_end = '2010-02-27T07:30:00.000' % string
%  ifplot = T/F; logical
%  
%
% Z. Eilon 04/2025
arguments
    nwk
    sta
    loc = '00';
    chan = 'BH*';
    time_start = '2010-03-27T00:00:00';
    time_end = '2010-03-28T05:00:00';
    ifplot = true;
end

quality = 'M'; % get spectrum from miniseed data banks
spd = 24*3600;


% put times in right format
% date_start_serial = convertTo(datetime(date_start),'datenum');
% date_end_serial = convertTo(datetime(date_end),'datenum');


%% Do the spectrum request
% build URL for request
mustang_url_preamble = 'https://service.earthscope.org/mustang/noise-psd/1/query?target=';
mustang_url = [mustang_url_preamble,...
    sprintf('%2s.%s.%s.%3s.%1s&',nwk,sta,loc,chan,quality),...
       sprintf('starttime=%s&',time_start),...
       sprintf('endtime=%s&',time_end),...
       'format=text'];
% e.g. https://service.earthscope.org/mustang/noise-psd/1/query?target=IU.ANMO.00.BH*.M&starttime=2013-10-22T00:00:00&endtime=2013-10-22T23:59:59&format=text

fprintf('%s\n',mustang_url)


% get MUSTANG infor
musti = urlread(mustang_url);

%% parse MUSTANG data
% this is a 1xN character array. must parse into lines
N = numel(musti);
 %find end of lines
eols = find(double(musti)==10);
Nlines = length(eols);

% get targets
Ntargets = textscan(musti(ilc(1,eols)),'# selected targets (%f) :'); Ntargets = Ntargets{1};
% make output structure
psdstruct = struct('nwk',[],'sta',[],'loc',[],'chan',[],'target',[],'timestart_covered',[],'timeend_covered',[]);
% put in metadata
for itarg = 1:Ntargets
    target_info_line = musti(ilc(1+itarg,eols));
    target_info = strsplit(target_info_line(11:end),'.');
    psdstruct(itarg).nwk = target_info{1};
    psdstruct(itarg).sta = target_info{2};
    psdstruct(itarg).loc = target_info{3};
    psdstruct(itarg).chan = target_info{4};
    psdstruct(itarg).target = target_info_line(11:end);
end
try
psdstruct(1).timestart_covered = musti(ilc(Ntargets+7,eols,10,inf));
psdstruct(1).timeend_covered = musti(ilc(Ntargets+8,eols,10,inf));
catch
    keyboard
end

% now start looping through targets, then looping through psds
iline = Ntargets + 9;
for itarg = 1:Ntargets
    iline = iline+1;
    % confirm target number
    if itarg ~= eval(musti(ilc(iline,eols,10,inf)))
        error('somehow wrong target')
    end
    iline = iline + 4;
    
    % if first target, put in start times etc.
    Npsd = eval(musti(ilc(iline,eols,12,inf)));
    % prepare psd time
    psdstruct(itarg).start_dts = datetime(ones(Npsd,3));
    psdstruct(itarg).end_dts = datetime(ones(Npsd,3));

    iline = iline-1; % rewind one
    for ipsd = 1:Npsd
        iline = iline+3;
        % get start and end times for this psd
        psdstruct(itarg).start_dts(ipsd) = datetime(musti(ilc(iline,eols,8,inf)));
        psdstruct(itarg).end_dts(ipsd) = datetime(musti(ilc(iline+1,eols,6,inf)));
        iline = iline+3;
        % now reading data depends on whether we know how many frequencies
        % there are. We will after the first time through. May be different
        % for each target, though.
        if ipsd == 1  
            ff = [];
            aa = [];
            iff = 1;
            while ~strcmp('#',musti(ilc(iline,eols,1,1)))
                lp = textscan(musti(ilc(iline,eols)),'%f,%f');
                ff(iff,1) = lp{1};
                aa(iff,1) = lp{2};
                % increment
                iff = iff+1;
                iline = iline+1;
                if iline==Nlines; break; end
            end
            Nfrq = length(ff);
            % insert as appropriate into this target's values
            psdstruct(itarg).frq = ff;            
            psdstruct(itarg).psd_matrix = nan(Nfrq,Npsd);
        else
            lp = textscan(musti(ilc(iline,eols,1):ilc(iline-1+Nfrq,eols,inf)),'%f,%f',Nfrq);
            ff = lp{1}(:);
            aa = lp{2}(:);
            iline = iline+Nfrq;
        end
        psdstruct(itarg).psd_matrix(:,ipsd) = aa;
    end

end

% times = psdstruct(1).start_dts + 0.5*(psdstruct(1).end_dts - psdstruct(1).start_dts);




%% plot 
if ifplot
    figh = figure(34); clf
    for itarg = 1:Ntargets
        subplot(Ntargets,1,itarg)
        plot(1./psdstruct(itarg).frq,psdstruct(itarg).psd_matrix,'-k','linewidth',1.5)
        set(gca,'ydir','normal','xscale','log','xlim',[0.11 180],'layer','top','fontsize',15)
        ylabel(gca,'Power [dB]','FontSize',18,'FontWeight','bold')
        title(psdstruct(itarg).target,'FontSize',18,'FontWeight','bold')
    end
    xlabel('Period (s)','FontSize',18,'FontWeight','bold')
end

end

%% SUBFUNCTIONS
function line_char_indices = ilc(iline,eols,ichar0,ichar1)
% function to get indices of all characters on line "iline", or,
% optionally, only the characters "ichars"
% where lines are separated by carriage return chars in positions "eols"
    if iline == 1
        line_char_indices = [1:eols(iline)-1];
    elseif iline > eols
        error('Not that many lines')
    else
        try
        line_char_indices = [eols(iline-1)+1:eols(iline)-1];
        catch
            keyboard
        end
    end
    if nargin == 3, ichar1 = ichar0; end
    if nargin < 3
        % do nothing, keep all chars
    else
        ichar0 = min([ichar0,length(line_char_indices)]);
        ichar1 = min([ichar1,length(line_char_indices)]);
        line_char_indices = line_char_indices(ichar0:ichar1);
    end

end