function [f,pdf,y,figh] = zeFetch_MUSTANGpdf(nwk,sta,loc,chan,date_start,date_end,ifplot,option,noiseprofile)
% [f,pdf,y,figh] = zeFetch_MUSTANGpdf(nwk,sta,loc,chan,date_start,date_end,ifplot,option,noiseprofile)
% 
% Function to get summary spectrum for a station for some date range. This
% will give an output obtained from aggregating individual hourly PSDs
% together to give some sort of summary statistic. Options are to get out a
% full PDF matrix (along with corresponding frequency and amplitude_dB
% bins) or spectra corresponding to mean/medium/mode/percentile values. 
% Note, will always request quality 'M', unless you go into code to change.
%
% INPUTS
%  nwk  = two-character code
%  sta  = four(or five) character code
%  loc = '*' % two char code
%  chan = '?H?'; % three char code
%  date_start = '2010-02-27' % string, date only
%  date_end = '2010-02-27' % string, date only
%  ifplot = T/F; logical
%  option = 'pdf' or 'spectra'; string
%  noiseprofile = 'mean,mode,median,min,max,37' % csv string list of types
%       of summary spectrum to get. Will recognise either numerical
%       percentiles (e.g. '37') or some special words, with some
%       equivalences, e.g. min==0%, max==100%, median==50%.
%  
%
% Z. Eilon 04/2025
arguments
    nwk
    sta
    loc = '*';
    chan = 'BHZ';
    date_start = '2010-02-27';
    date_end = '2010-03-28';
    ifplot = true;
    option = 'spectra';
    noiseprofile = 'mean,min,mode,max';
end

quality = 'M'; % get spectrum from miniseed data banks
spd = 24*3600;


% put times in right format
date_start_serial = convertTo(datetime(date_start),'datenum');
date_end_serial = convertTo(datetime(date_end),'datenum');


%% Do the spectrum request
% build URL for request
mustang_url_preamble = 'http://service.earthscope.org/mustang/noise-pdf/1/query?target=';
mustang_url_base = [mustang_url_preamble,...
    sprintf('%2s.%s.%s.%3s.%1s&',nwk,sta,loc,chan,quality),...
       sprintf('start=%s&',date_start),...
       sprintf('end=%s',date_end),...
        ];
% e.g. https://service.earthscope.org/mustang/noise-pdf/1/query?target=IU.ANMO.00.BHZ.M&starttime=2013-12-01&endtime=2013-12-31&format=text
% e.g. https://service.earthscope.org/mustang/noise-pdf/1/query?target=IU.ANMO.00.BHZ.M&starttime=2010-12-01&endtime=2017-12-31&noiseprofile.type=min,mean,mode,median,max&format=noiseprofile_text

% now add on format things per the inputs
switch option
    case 'pdf'
        mustang_url = [mustang_url_base,'&format=text'];
    case 'spectra'
        mustang_url = [mustang_url_base,'&noiseprofile.type=',noiseprofile,'&format=noiseprofile_text'];        
        nspec = sum(noiseprofile==',')+1;
end
fprintf('%s\n',mustang_url)


% get MUSTANG data
mustangdat = urlread(mustang_url);

% parse MUSTANG data
[f,pdf,y] = parse_mustangdat_pdf(mustangdat);

if ishandle(ifplot)
    figh = ifplot; 
    figure(figh);
    ifplot = true;
elseif ifplot
    figh = figure;
end

if ifplot
    switch option
        case 'pdf'
            imagesc(1./f,y,pdf);
            grid on
            colormap(jet)

            % put on noise models
            hold on
            okfnm = (f < 10) & (f > 0.001);
            plot(1./f(okfnm),nhnm(f(okfnm)),'-','linewidth',2,'Color',0.7*[1 1 1])
            plot(1./f(okfnm),nlnm(f(okfnm)),'-','linewidth',2,'Color',0.7*[1 1 1])
            
        case 'spectra'
            hp = plot(1./f,pdf,'linewidth',2.5);
            legend(hp,string(y),'fontweight','bold');
    end

    set(gca,'ydir','normal','xscale','log','xlim',[0.05 200],'layer','top','fontsize',15)
    xlabel('Period (s)','FontSize',18,'FontWeight','bold')
    ylabel('Power [dB]','FontSize',18,'FontWeight','bold')

            

end

end

%% SUBFUNCTIONS
function [f,pdf,values,starttime,endtime] = parse_mustangdat_pdf(mustangdat)
% this is a 1xN character array. must parse into lines
N = numel(mustangdat);
 %find end of lines
eols = find(double(mustangdat)==10);
Nlines = length(eols);
Nf = Nlines-5;
% get indices of characters corresponding to each block of information
starttime_ichars = [eols(1)+2:eols(2)-1]; % eol of previous line, + 1 to get to current line, +1 to ignore '#' character
endtime_ichars = [eols(2)+2:eols(3)-1];   % eol of previous line, + 1 to get to current line, +1 to ignore '#' character 
header_ichars = [eols(4)+2:eols(5)-1];    % eol of previous line, + 1 to get to current line, +1 to ignore '#' character 
data_ichars = [eols(5)+1:eols(end)-1]; % assume last char is an eol.

% extract metadata
starttime = datetime(regexprep(mustangdat(starttime_ichars),'start=',''));
endtime = datetime(regexprep(mustangdat(endtime_ichars),'end=',''));
header = strsplit(mustangdat(header_ichars),', ');
values = header(2:end);
% turn percentiles into numbers from strings
for iv = 1:length(values)
    if any(values{iv}=='%')
        values{iv} = eval(regexprep(values{iv},'%',''));
    end
end
Nval = length(header);

% extract data
dataformat = repmat('%f,',1,Nval);
data = textscan(mustangdat(data_ichars),dataformat(1:end-1));


% get out frequencies (either with all repeats for pdf, or unique.
f = data{1};
% get out spectral data
pdf = nan(Nf,Nval-1);
for iv = 1:Nval-1
    pdf(:,iv) = data{iv+1};
end

if strcmp( values{end},'hits') % we are looking at a PDF
    ff = unique(f); 
    a = pdf(:,1);
    hits = pdf(:,2);
    aa = unique(a); % old way - only unique amplitude bins
    aa = [min(aa):max(aa)]'; % new way - all amplitudes in range
    Nf = length(ff);
    Na = length(aa);
    % Note, file only has entries for non-zero values - have to loop
    pdf = zeros(Na,Nf);
    for iff = 1:Nf
        thisf = find(f==ff(iff));
        [~,IA] = intersect(aa,a(thisf));
        pdf(IA,iff) = hits(thisf);
    end
    % and now organize outputs
    values = aa;
    f = ff;
else % we are looking at spectra - no need to alter

end





end