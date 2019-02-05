function plot_resp( resp,faxis,extra_zp )
% plot_resp( resp,faxis,extra_zp )
% function to plot instrument response from data given poles, zeros and gain
% adapted by Zach Eilon from function written by Ge Jin, 2014/02/27
% 
% INPUTS 
%   resp     = complex response vector
%   zeros    = frequency vector
%   extra_zp = [0] flag to say that there are extra zeros (positive) or poles (negative) 
%              in the response, which should be accounted for when plotting the response
if nargin<3
    extra_zp = 0;
end

w = 2*pi*faxis;
resp_plot = resp./(w.^extra_zp);

figure(35), clf, set(gcf,'position',[360   514   900   800]);
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

