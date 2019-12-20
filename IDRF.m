function [RF_Time, RF] = IDRF(P,D,dt,t_bounds,gauss_t,accept_mis,itmax)
% [RF_Time, RF] = IDRF(P,D,dt,t_bounds,gauss_t,accept_mis,itmax)
%
% Iterative Deconvolution and Receiver-Function Estimation in time domain
%
% P for parent phase (S in Sp case), D for daughter phase.
% dt = time delta (1/samprate) (s)
% t_bounds = time bounds to preserve in eventual RF (really just cut out past t_max)
% gauss_t = 1 std for the gaussian convolved with the delta spikes (recommend 1 s)
% accept_mis = acceptable misfit for (D - conv(P,rf))/length(D)
% itmax = max number of iterations
%
% Sent by Karen Fischer, Jan 2019 
% I think developed by Nick Mancinelli in 2017
% adapted by Hannah Krueger and/or Junlin Hua, early 2019
% lightly edited by Z. Eilon, Jan 2019

t_max = max(t_bounds); %max preserved point in RF\
t_min = min(t_bounds); %min preserved point in RF\

misfit = 999999999999999999;
% misfit = sqrt(sum((P-D).^2)); << was

RF = zeros(length(P)*2-1,1);

D_cur = D;

itnum = 0;


% while (itnum <= itmax) && (misfit_old > accept_mis) << was
while (itnum <= itmax) && (misfit > accept_mis)
    [amp_corr,t_corr] = xcorr(D_cur,P);
    auto_corr = xcorr(P);
    [~,ind] = max(abs(amp_corr));
    amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
    RF(ind) = RF(ind)+amp_rf;
    D_sub = conv(P,RF,'same');
    D_cur = D - D_sub;
    %plot(D_cur)
%     misfit_old = misfit; << was
    misfit = sqrt(sum(D_cur.^2))/length(D_cur);
    itnum = itnum+1;
end



RF_Time = t_corr*dt;

RF(RF_Time>t_max)=0;
RF = RF(RF_Time<=t_max);
RF_Time = RF_Time(RF_Time<=t_max);
% RF(RF_Time>t_max & RF_Time<t_min)=0;
% RF = RF(RF_Time<=t_max & RF_Time>=t_min);
% RF_Time = RF_Time(RF_Time<=t_max & RF_Time>=t_min);

if gauss_t~=0
    % gauss_len = length(RF);
    % Gauss_win = gausswin(gauss_len,gauss_sig*4);
    gauss_sig = gauss_t/dt;
    x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
    Gauss_win = exp(-x.^2/(2*gauss_sig^2));
    RF = conv(RF,Gauss_win,'same');
end

RF = flipud(RF);
RF_Time = fliplr(RF_Time);

%plot(dt*t_corr(t_corr<0 & dt*t_corr>-40),RF_a(t_corr<0 & dt*t_corr>-40))




end