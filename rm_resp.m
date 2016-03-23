function odat = rm_resp(idat,zz,pp,gain,samprate,extra_zp,ifplot)
%% function to remove instrument response from data given poles, zeros and gain
% adapted by Zach Eilon 2015/09/10 from function written by Ge Jin, 2014/02/27
% 
% INPUTS 
%   idat     = data matrix, with time series in columns [ N x ntraces ] 
%   zz       = vector of complex zeros
%   pp       = vector of complex poles
%   gain     = gain amplitude
%   samprate = sample rate, in Hz
%   extra_zp = [0] flag to say that there are extra zeros (positive) or poles (negative) 
%              in the response, which should be accounted for when plotting the response
%   ifplot   = [0] option to plot the response function,
%              and the data before & after. 
%
% OUTPUTS
%   odat     = corrected data matrix, with time series in columns [ N x ntraces ] 

if nargin<7
    ifplot=0;
end
if nargin<6
    extra_zp=0;
end

lo_corner = 0.005;  % in Hz - was 0.005
npoles_myfilt=5; 

N = size(idat,1);
delta = 1/samprate;
T = N*delta;

data = detrend(idat);
data = flat_hanning_win(1:N,data,1,N,50);

% make freq. axis
if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

% calc response function in f-domain
w = faxis.*2*pi;
resp = ones(size(w));
for ip = 1:length(pp)
	resp = resp./(1i*w - pp(ip));
end
for iz = 1:length(zz)
	resp = resp.*(1i*w - zz(iz));
end
resp = resp*gain;

resp_plot = resp./(w.^extra_zp);

% add extra high-pass filter to get rid of >200s stuff - Jingle's input
lo_w=2*pi*lo_corner;
hpfiltfrq=( ((w./lo_w).^(2*npoles_myfilt))./(1+(w./lo_w).^(2*npoles_myfilt)) );
norm_trans=hpfiltfrq./resp;    % this is normalization transfer function
norm_trans(isnan(norm_trans)) = 0;

% Option to put on a low-pass filter for the high end stuff?! Prob not
% needed.

fftdata = fft(data);
fftdata = fftdata(:).*norm_trans(:);
odat = real(ifft(fftdata));

if ifplot
figure(39), clf, set(gcf,'position',[360   514   900   800]);

ind = faxis>0;

% amplitude response
subplot(3,2,1), hold on, set(gca,'fontsize',18,'xscale','log','yscale','log','FontSize',14)
loglog(faxis(ind),abs(resp_plot(ind)),'r','Linewidth',1.5);
xlim([1/T 2/delta])
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Response','FontSize',18)

% phase response
subplot(3,2,2), hold on, set(gca,'fontsize',18,'xscale','log','Ytick',[-180:90:180],'FontSize',14)
plot(faxis(ind),r2d(angle(resp_plot(ind))),'r','Linewidth',1.5);
xlim([1/T 2/delta]), ylim([-180 180])
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Response','FontSize',18)

% raw data
subplot(3,2,3:4), hold on, set(gca,'fontsize',14)
plot(delta*(0:N-1),idat,'b');

% response-removed data
subplot(3,2,5:6), hold on, set(gca,'fontsize',14)
plot(delta*(0:N-1),odat,'g');
xlabel('Time (s)','FontSize',18)

end

return

