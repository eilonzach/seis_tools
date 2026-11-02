function [tt,aa,info] = zeFetch_MUSTANGtimeseries(nwk,sta,loc,chan,date_start,date_end,periods,ifavg,ifplot)
% [tt,aa] = zeFetch_MUSTANGtimeseries(nwk,sta,loc,chan,date_start,date_end,periods,ifavg,ifplot)
% 
% Function to get time series of psd values for given frequency range. Can
% only handle a single "target", so must do one nwk+sta+loc+chan at a time.
% Note, will always request quality 'M', unless you go into code to change.
% Can put in a range of periods, and can also put in two periods and then
% take the average of all the periods in between those with ifavg.
%
% INPUTS
%  nwk  = two-character code
%  sta  = four(or five) character code
%  loc = '*' % two char code
%  chan = '?H?'; % three char code
%  time_start = '2010-02-27T06:30:00.000' % string
%  time_end = '2010-02-27T07:30:00.000' % string
%  periods = [minperiod,maxperiod] % or a list of several periods. double
%  ifavg = T/F; logical, option to average over all the outputs
%  ifplot = T/F; logical  
%
% Z. Eilon 04/2025
arguments
    nwk
    sta
    loc = '00';
    chan = 'BHZ';
    date_start = '2010-01-01';
    date_end = '2014-01-01';
    periods = [5 9];
    ifavg = true;
    ifplot = true;
end

quality = 'M'; % get spectrum from miniseed data banks
spd = 24*3600;


% put times in right format
% date_start_serial = convertTo(datetime(date_start),'datenum');
% date_end_serial = convertTo(datetime(date_end),'datenum');


%% Do the spectrum request
% information about the periods
if ~isempty(periods)
    pstr = '&periods=[';
    for ip = 1:length(periods)
        pstr = [pstr,num2str(periods(ip))];
        if ip~=length(periods)
            pstr = [pstr,','];
        end
    end
    pstr = [pstr,']'];
else 
    pstr = '';
end

% build URL for request
mustang_url_preamble = 'https://service.earthscope.org/mustang/noise-mode-timeseries/1/query?target=';
mustang_url = [mustang_url_preamble,...
    sprintf('%2s.%s.%s.%3s.%1s&',nwk,sta,loc,chan,quality),...
       sprintf('starttime=%s&',date_start),...
       sprintf('endtime=%s&',date_end),...
       'format=text',pstr];
% e.g. https://service.earthscope.org/mustang/noise-psd/1/query?target=IU.ANMO.00.BH*.M&starttime=2013-10-22T00:00:00&endtime=2013-10-22T23:59:59&format=text

fprintf('%s\n',mustang_url)

% get MUSTANG info
musti = urlread(mustang_url);

%% parse MUSTANG data
% this is a 1xN character array. must parse into lines
N = numel(musti);
 %find end of lines
eols = find(double(musti)==10);
Nlines = length(eols);

% make + fill output structure
info = struct('nwk',[],'sta',[],'loc',[],'chan',[],'target',[],...
    'timestart_requested',[],'timeend_requested',[],...
    'timestart_returned',[],'timeend_returned',[]);

target_info = strsplit(musti(ilc(1,eols,9,inf)),'.');
request_info = strsplit(musti(ilc(2,eols,23,inf)),',');
return_info = strsplit(musti(ilc(3,eols,22,inf)),',');
info.nwk = target_info{1};
info.sta = target_info{2};
info.loc = target_info{3};
info.chan = target_info{4};
info.target = musti(ilc(1,eols,9,inf));
info.timestart_requested = datetime(request_info{1});
info.timeend_requested = datetime(request_info{2});
info.timestart_returned = datetime(return_info{1});
info.timeend_returned = datetime(return_info{2});

% get N periods
Np = eval(musti(ilc(5,eols,8,inf)));
% work out N days
Nday = (Nlines-3)./(3+Np);

% prep output
tt = datetime(zeros(Nday,3));
aa = nan(Nday,Np);
% ff = nan(Np);
info.pp = nan(Np);

% now start looping through days
iline = 4;
for id = 1:Nday
    tt(id) = datetime(musti(ilc(iline,eols,6,inf)));
    iline = iline + 3;
    lp = textscan(musti(ilc(iline,eols,1):ilc(iline-1+Np,eols,inf)),'%f,%f',Np);
    ff = lp{1}(:);
    aa(id,:) = lp{2}(:);
    iline = iline+Np;
    if id==1
        info.pp = 1./ff;
    end
end

% now do averaging if desired
if ifavg
    aa = mean(aa,2);
    info.pp = mean(info.pp);
end



%% plot 
if ifplot
    figh = figure(35); clf
    plot(tt,aa,'-k','linewidth',1.5)
    hold on 
    plot(tt,smooth(aa,30),'-m','linewidth',2)
    set(gca,'ydir','normal','layer','top','fontsize',15)
    ylabel(gca,'Power [dB]','FontSize',18,'FontWeight','bold')
    xlabel('Time','FontSize',18,'FontWeight','bold')
    title(info.target,'FontSize',18,'FontWeight','bold')
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