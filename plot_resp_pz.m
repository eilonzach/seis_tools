function plot_resp_pz( pp,zz,gain, fmin, fmax,extra_zp, fign )
% plot_resp_pz( poles,zer0s,gain,[fmin],[fmax],[extra_zp])
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

if nargin<4  || isempty(fmin)
    fmin = 1e-4;
end
if nargin<5 || isempty(fmax)
    fmax = 1e3;
end
if nargin<6
    extra_zp = 0;
end
if nargin<7
    fign=33;
end

if extra_zp~=0
    fprintf('Removing %.0f zeros from response\n',extra_zp)
    if extra_zp > 0
        zz = sort(zz);
        zz = zz(1+extra_zp:end);
    else
        zz = [zeros(abs(extra_zp),1);zz];
    end
end


fprintf('%-2.0f Zeros (in log10 f-space) at:\n',length(zz))
for jj = 1:length(zz)
    fprintf('Zero %2.0f:  %9.3e %10.3ei  (freq. = %.3f Hz)\n',jj,real(zz(jj)),imag(zz(jj)),abs(zz(jj)))
end
% log10(abs(zz)./(2*pi))
fprintf('%-2.0f Poles (in log10 f-space) at:\n',length(pp))
for kk = 1:length(pp)
    fprintf('Pole %2.0f:  %9.3e %10.3ei  (freq. = %.3f Hz)\n',kk,real(pp(kk)),imag(pp(kk)),abs(pp(kk)))
end
% log10(abs(pp)./(2*pi))

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

resp_plot = resp;%./(1i*w.^extra_zp);

figure(fign), clf, set(gcf,'position',[360   514   900   800],'color','w');
% amplitude 
subplot(211), hold on, 
loglog(faxis,abs(resp_plot),'r','Linewidth',2);
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Response amplitude','FontSize',18)
set(gca,'fontsize',18,'xscale','log','yscale','log','FontSize',14,'box','on','linewidth',1.5), grid on

% phase 
subplot(212), hold on, 
plot(faxis,r2d(angle(resp_plot)),'r','Linewidth',2);
ylim([-180 180])
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Response phase','FontSize',18)
set(gca,'fontsize',18,'xscale','log','Ytick',[-180:90:180],'FontSize',14,'box','on','linewidth',1.5), grid on

end

