function [ datwf,datf,datwc,datc,fdeets,ttws,tts ] = data_clean( traces,cleaning_parm,model )
% [ datwf,datf,datwc,datc,fdeets,ttsw,tts ] = data_clean( traces,cleaning_parm )
%Function to clean (detrend etc.) filter and window data from N stations
% Taper and window chosen so that prex=0, postx=10, taperx = 0.1 will be
% start at prex (=0), ramp up to 1 at nwin/taperx (=1), and then ramp down from 1 to 0
% between postx - nwin/taperx (=9) and posttime (=10)
%
% INPUTS
%  traces = [npt x nsta] matrix of the data in columns
%  cleaning_parameters = structure with parameters for windowing and
%   filtering etc. with fields:
%     cleaning_parm.samprate = sample rate per second
%     cleaning_parm.pretime  = time before 0s that traces start
%     cleaning_parm.prex     = time before 0s for window to start
%     cleaning_parm.postx    = time after 0s for window to end
%     cleaning_parm.taperx   = proportion of trace that is taper (each end)
%     cleaning_parm.fhi      = high bp freq.
%     cleaning_parm.flo      = low bp freq.
%     cleaning_parm.filtopt  = low bp freq.
%     cleaning_parm.npoles   = number of poles for filter
%     cleaning_parm.npass    = # of passes for filter: 1=causal, 2=acausal
%     cleaning_parm.norm     = whether (1) or not (0) to normalise traces by std
%     cleaning_parm.detrend  = whether (1) or not (0=default) to detrend traces
%
% OUTPUTS
%  datwf   = cleaned, windowed, tapered, filtered data in columns
%  datf    = cleaned, filtered data in columns
%  datwc    = cleaned, windowed, tapered data in columns
%  datc    = cleaned data in columns
%  filtnam = string with filter information
%  ttws    = vector of times for the filtered, windowed traces
%  tts     = vector of times for unfiltered, unwindowed, clean traces
%
% Written by Z. Eilon 08/2015

cp = cleaning_parm;

nsta = size(traces,2);
npt = size(traces,1);
dt = 1./cp.samprate;

if ~isfield(cp,'npass'), cp.npass = 2; end


%% Build filters
fmin = cp.samprate/npt;
fmax = 0.5/dt;
% edit filter frequencies in case they don't work...
% flo = max([cp.flo,fmin]); %max period can't be longer than window length
% fhi = min([cp.fhi,fmax]); %min period can't be shorter than Nyquist period
flo = cp.flo;
fhi = cp.fhi;

fdeets = [fhi flo cp.npoles];
% ffund = 1./(npt.*dt);
% % option 1: butter
% if flo > ffund && fhi.*dt.*2<1
%     [bb,aa]=butter(cp.npoles, [flo, fhi].*dt.*2);
% elseif flo > ffund && fhi.*dt.*2==1
%     [bb,aa]=butter(cp.npoles, flo.*dt.*2,'high');
% elseif  flo<=ffund && fhi.*dt.*2<1
%     [bb,aa]=butter(cp.npoles, fhi.*dt.*2,'low');
% elseif  flo<=ffund && fhi.*dt.*2==1
%     datf = traces;
%     return
% end
% % option 2: cheby
% [bb,aa]=cheby1(cp.npoles,0.5, [flo, fhi].*dt.*2);
% % option 3: higher order butter
% [z,p,k]=butter(cp.npoles, [flo, fhi].*dt.*2.);
% [sos,g]=zp2sos(z,p,k); bb=sos; aa=g;

if ~isfield(cp,'detrend')
    cp.detrend = 0;
end

%% Make the taper window
% WAS nwin=round((cp.postx+cp.prex)/dt)+1; % window length in samples
nwin=round((cp.postx+cp.prex)/dt); % window length in samples
n1=round((cp.pretime-cp.prex)/dt); % first sample in window
wdo2=[zeros(n1,1);tukeywin(nwin,2*cp.taperx); zeros(npt-n1-nwin,1)]; % taperx% tukey window
ibds=[max(1,n1-floor(nwin*2*cp.taperx)), min(npt,n1+nwin+floor(nwin*2*cp.taperx))]; % extend so traces are extra 20% of window on either side
jbds=ibds(1):ibds(2); % indices of points to keep

datwf=zeros((ibds(2)-ibds(1)+1), nsta);
datwc=zeros((ibds(2)-ibds(1)+1), nsta);
datf=zeros(size(traces));
datc=zeros(size(traces));

for is=1:nsta
    %% Clean traces
    rec=traces(:,is);
    % isfi=find(isfinite(rec));
    % mn=nanmean(rec(isfi));    % Deal with NaNs
    % inan=find(isnan(rec)); % find NaNs
    nnan=find(~isnan(rec)); % find not-NaNs
    if isempty(nnan), datwf(:,is) = NaN; continue, end
    % rec(nnan)=rec(nnan)-mn; % take off mean 
    if cp.detrend
        rec(nnan)=detrend(rec(nnan)); % detrend not-nans
    end
    rec(isnan(rec))=0; %set NaNs to zero
    wdo1 = [zeros(nnan(1)-1,1);tukeywin(length(nnan),2*cp.taperx);zeros(length(rec)-nnan(end),1)];
    rec = rec.*wdo1; % first taper of whole window. Sets NaNs to zero and tapers end of non-NaN
    if (length(nnan)<2)      % skip over garbage trace
      fprintf('Garbage trace number %d\n');
      datc(:,is)=zeros(size(rec));
      continue
    end
    datc(:,is)=rec;
    
    % window
    if nargout>2
    try
    recw=rec.*wdo2;
    datwc(:,is)=recw(jbds);
    catch
        fprintf('ERROR WITH WINDOWING IN DATA CLEAN\n')
        save('traces','traces'),save('cp','cp'),try save('model','model');end
    end
    end

    %% Filtered data
%     % pad with plenty of zeros for the filter
%     rec = [zeros(1000,1);rec;zeros(1000,1)]; 
%         % option 1: filter with phase
%         % recf=filter(bb, aa, rec);
%     % zerophase bandpass filter & window
%     recf=filtfilt(bb, aa, rec);
%     % lop off padding
%     recf = recf(1001:end-1000);
    recf = filt_quick(rec,flo,fhi,dt,cp.npoles,cp.npass);

    % detrend again
    if cp.detrend
        recf=detrend(recf);
    end

    datf(:,is)=recf;
    

    % window
    try
    recwf=recf.*wdo2;
    datwf(:,is)=recwf(jbds);
    catch
        fprintf('ERROR WITH WINDOWING IN DATA CLEAN\n')
        save('traces','traces'),save('cp','cp'),try save('model','model');end
    end
    
    % normalize
    if cp.norm==1
    	datf(:,is)=datf(:,is)./std(recf);
        datwf(:,is)=datwf(:,is)./std(recf);
    end

end % stas loop

tts =  (0:(npt-1)).*dt - cp.pretime;        tts = tts(:);  % columnise
ttws = ([ibds(1):ibds(2)]-1).*dt - cp.pretime;  ttws = ttws(:); % columnise

end % on function
