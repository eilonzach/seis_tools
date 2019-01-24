function [RF_Time,RF] = IDRF(P,D,dt,t_cutoff,t_trun,gauss_t,min_misf_fchange,min_misf_fdrop,nitermax)
% [RF_Time,RF] = IDRF(P,D,dt,t_for,t_trun,gauss_t,accept_mis1,accept_mis2,itmax)
% Iterative Deconvolution and Receiver-Function Estimation in time domain
% P for parent phase (S in Sp case), D for daughter phase.
% gauss_t: 1 std for the Gaussian convolved
% accept_mis1, accept_mis2 and itmax for stop criteria of the loop
% t_for,t_trun: remove amplitude with t > 0 in the Sp case 
% 
% Sent by Karen Fischer, Jan 2019 
% I think developed by Nick Mancinelli in 2017
% adapted by Z. Eilon, Jan 2019

misfit_old = 99999999;
misfit_ref = sqrt(sum(D.^2));
misfit = misfit_ref;

%RF = zeros(length(P)*2-1,1);
RF_tmp = zeros(length(P)*2-1,1);

D_cur = D;

itnum = 0;

[~,t_corr] = xcorr(D_cur,P);
RF_Time = t_corr*dt;

while misfit_old-misfit>min_misf_fchange*misfit && misfit>min_misf_fdrop*misfit_ref && itnum <= nitermax
    [amp_corr,~] = xcorr(D_cur,P);
    auto_corr = xcorr(P);
    [~,ind] = max(abs(amp_corr));
    amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
    %
    %                 if RF_Time(ind) > -9.5 && (amp_rf>0 || abs(amp_rf)<0.02)
    %                     RF_tmp(ind) = RF_tmp(ind)+amp_rf;
    %                 else
    %                     RF_tmp(ind) = RF_tmp(ind)+amp_rf;
    %                     RF(ind) = RF(ind)+amp_rf;
    %                 end
    RF_tmp(ind) = RF_tmp(ind)+amp_rf;
    D_sub = conv(P,RF_tmp,'same');
    D_cur = D - D_sub;
    %plot(D_cur)
    misfit_old = misfit;
    misfit = sqrt(sum(D_cur.^2))/misfit_ref;
    itnum = itnum+1;
end

RF = RF_tmp;


RF(RF_Time>t_trun)=0;
RF = RF(RF_Time<=t_cutoff);
RF_Time = RF_Time(RF_Time<=t_cutoff);

if gauss_t~=0
    gauss_sig = gauss_t/dt;
    x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
    Gauss_win = exp(-x.^2/(2*gauss_sig^2));
    RF = conv(RF,Gauss_win,'same');
    
end
%RF = flipud(RF);
%RF_Time = fliplr(RF_Time);
end