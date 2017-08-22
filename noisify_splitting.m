function [ odat,SNRest,dofps ] = noisify( idat,ndat,samprate,SNR,plotopt,theta,phi,dT )
% function [ odat,SNRest,dofps ] = noisify( idat,ndat,samprate,SNR,plotopt,theta,phi,dT )
% 
% Function to take an input set of data (idat) and add noise with the power
% distribution of an input noise series (ndat). Ideally the two input
% series should be of the same length and sample rate
%
% INPUTS:
% idat     - time series of data in columns (e.g. 3 channels would be an
%              nsamps x 3 matrix) (try to have E,N,Z)
% ndat     - time series of noisy data in columns (e.g. 3 channels would be
%              an nsamps x 3 matrix). Noise window should have no arrival. 
% samprate - sample rate of input time series (samples per second)
% SNR      - signal to noise ratio desired
% plotopt  - 1 to plot, anything else will not.
% -- if you have them --
% theta    - forward azimuth of incoming SKS wave (relative to north)
% phitrue  - true fast direction (relative to north)
% dttrue   -  true time delay between fast and slow
% 
% OUTPUTS:
% odat     - output time series, same size as idat and with noise added
% SNRest   - measured SNR of odat
% dofps    - degrees of freedom per second of the output series
global SNRest dSNR

% disp('start')
if nargin < 5
    plotopt=0;
end

%% Key values
noise=1./SNR;
dt=1./samprate;
nsamps=size(idat,1); 
nchans=size(idat,2);
samplength=nsamps/samprate;
if size(ndat,1)~=nsamps || size(ndat,2)~=nchans 
    error('Make idat and ndat the same size!')
end

%% Detrend and scale input time series
ndat=detrend(ndat,'constant'); %detrend
idat=detrend(idat,'constant'); %detrend
%%%%%%%%% adjust amplitudes %%%%%%%%%%%%%%
for ic=1:nchans; ndat(:,ic)=ndat(:,ic)./sqrt(mean(ndat(:,ic).^2)); end 
% for ic=1:nchans; idat(:,ic)=idat(:,ic)./sqrt(mean(idat(:,ic).^2)); end 


%% GET PSD of NOISE
% Window and get sample of noise
% get psd for each channel
if exist('psde','var')==1
    disp('N.B. Using existing psd')
else
nw=4;
psde=zeros(nsamps/2 + 1,4);
freq=[1:floor(nsamps/2)+1]'./(nsamps*dt);%in mHz - to get in Hz, divide by 1000
for ic=1:nchans
[Pxx,~] = pmtm(ndat(1:nsamps,ic),nw,freq,samprate);
H = dspdata.psd(Pxx,'Fs',samprate);
if plotopt==1
figure(10)
plot(H)
end
psde(:,1)=H.Frequencies;
psde(:,ic+1)=H.Data;
end
end

%% Put in noise 
% make fft of synth traces and add to fft of a bootstrap noise series
%%%% REPEAT OVER NOISE REALISATIONS UNTIL SNRest IS CLOSE TO SNR IN %%%%
dSNR=1000; 
while dSNR > SNR/10 
%     disp('go')
clear('i');
FA=zeros(nsamps,3);
FB=zeros(nsamps,3);
% Ft of noise
for ic=1:nchans
FA(1,ic)=psde(1,ic+1);
FA(nsamps/2+1,ic)=psde(end,ic+1);
FA(2:nsamps/2,ic)=psde(2:end-1,ic+1).*exp(1i*2*pi*random('unif',0,1,nsamps/2-1,1));
FA(nsamps/2+2:end,ic)=flipud(conj(FA(2:nsamps/2,ic)));
FA(:,ic)=FA(:,ic)/sqrt(mean(abs(FA(:,ic)).^2));
FA(:,ic)=FA(:,ic)*sqrt(nsamps);
end
% Ft of signal
for ic=1:nchans
FB(:,ic)=fft(idat(:,ic));
end
% Add the two Fts, scaling noise to SNR
FF=(noise*FA + FB);
% FF=(FB);

%ifft to get data:
odat=zeros(nsamps,3);
for ic=1:nchans
    odat(:,ic)=ifft(FF(:,ic));
    if plotopt==1
    figure(13);
	A = max(max(idat));
    subplot(nchans,2,2*ic); plot([0:dt:samplength-dt],odat(:,ic));
    ylim([-A  A])
    subplot(nchans,2,2*ic-1); plot([0:dt:samplength-dt],idat(:,ic));
    ylim([-A  A])
    ylabel(sprintf('Component %u',ic));
    end
end

%% IF POSS, UNSPLIT AND MEASURE SIGNAL TO NOISE
if nargin > 5
    
dtheta = phi-theta; %I think this is right - this is the clockwise angle by which the true radial is rotated to get to the fast
shiftt = ceil(dT/dt);
dphi = -phi;

de=odat(:,1);
dn=odat(:,2);
%taper
  len  = round(length(de)*.05); %taper length is 3% of total seismogram length
  nn   = 1:len;
  nn2  = (length(de)-len+1):length(de);
  x    = linspace(pi, 2*pi, len);
taper  = 0.5 * (cos(x')+1);
taper2 = fliplr(taper);
de(nn) = de(nn).*taper;     de(nn2) = de(nn2).*taper2;
dn(nn) = dn(nn).*taper;     dn(nn2) = dn(nn2).*taper2;
% filter NB THIS MIGHT NEED TO BE CHANGED DEPENDING ON THE NOISE 
[b,a]  = butter(3, 2*[0.02 0.125]/samprate);
xnf=filtfilt(b,a,dn');
xef=filtfilt(b,a,de');
% unsplit
cr=cosd(-dphi);
sr=sind(-dphi);
rrmat=[cr,sr;-sr,cr];
xf=rrmat(1,:)*[xnf;xef]; xf=xf'; 
xs=rrmat(2,:)*[xnf;xef]; xs=xs';
xf=xf(1:end-shiftt);
xs=xs(shiftt+1:end);

cr=cosd(-dtheta);
sr=sind(-dtheta);
rrmat=[cr,sr;-sr,cr];
xq=rrmat(1,:)*[xf';xs']; xq=xq'; 
xt=rrmat(2,:)*[xf';xs']; xt=xt';
SNRest=max(abs(xq)) / (2*std(xt));
dSNR=abs(SNRest-SNR);
else % if no unsplit parameters available
    dSNR=0;
    for ic=1:nchans
        SNRest(ic)=max(abs(odat(:,ic)))/sqrt(mean(odat(:,ic).^2));
    end
end % if nargin includes unsplit parameters

end % while dSNR too big

%% Degrees of freedom per second
% calculate using first zero-crossing of odat acf
dofps = zeros(nchans,1);
for ic = 1:nchans
    spdof = min(abs(zerof(xcorr(odat(:,ic))) - nsamps))/samprate;
    dofps(ic) = 1./spdof;
end

