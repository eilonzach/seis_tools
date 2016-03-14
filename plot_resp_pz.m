function plot_resp_pz( pp,zz,gain, fmin, fmax,extra_zp, fign )
%plot_resp_pz( poles,zer0s,gain,[fmin],[fmax],[extra_zp])
% function to plot instrument response from data given poles, zeros and gain
% written by Z. Eilon -- 09/2015
% 
% INPUTS 
%   pp       = vector of complex poles
%   zz       = vector of complex zeros
%   gain     = gain amplitude
%   fmin     = [1e-4] minimum freq. (i.e. 1/T where T is windlength in s)
%   fmax     = [1e2] maximum freq. (i.e. 2/dt, Nyquist freq.) 
%   eztra_zp = [0] flag to say that there are extra zeros (positive) or poles (negative) 
%              in the response, which should be accounted for when plotting the response
%   fign     = [33] figure number

if nargin<4
    fmin = 1e-4;
end
if nargin<5
    fmax = 1e3;
end
if nargin<6
    extra_zp = 0;
end
if nargin<7
    fign=33;
end

fprintf('Poles (in log10 f-space) at:\n')
log10(abs(pp)./(2*pi))
fprintf('Zeros (in log10 f-space) at:\n')
log10(abs(zz)./(2*pi))

% make freq. axis
faxis = [0,fmin:fmin:fmax]';

% calc response function in f-domain
w = faxis.*2*pi;

%% Calc. response
resp = ones(size(w));
for ip = 1:length(pp)
	resp = resp./(1i*w - pp(ip));
end
for iz = 1:length(zz)
	resp = resp.*(1i*w - zz(iz));
end
resp = resp*gain;

resp_plot = resp./(1i*w.^extra_zp);

figure(fign), set(gcf,'position',[360   514   900   800]);
% amplitude 
subplot(211), hold on, set(gca,'fontsize',18,'xscale','log','yscale','log','FontSize',14), grid on
loglog(faxis,abs(resp_plot),'r','Linewidth',1.5);
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Response amplitude','FontSize',18)

% phase 
subplot(212), hold on, set(gca,'fontsize',18,'xscale','log','Ytick',[-180:90:180],'FontSize',14), grid on
plot(faxis,r2d(angle(resp_plot)),'r','Linewidth',1.5);
ylim([-180 180])
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Response phase','FontSize',18)


end

