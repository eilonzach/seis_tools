function [ dat_down ] = downsamp( dat, samprate, resamprate )
%[ dat_downsamp ] = downsamp( dat, samprate, resamprate )
%   Function to downsample a data series from samprate to resamprate. If
%   necessary, apply a low-pass filter with the resamprate as the nyquist
%   in order to prevent aliasing. dat can be a column vector or a matrix
%   with data in columns.
%   0.5% of time window at each end will be tapered. Deal with it. 

nsamps = size(dat,1);

tt = [0:nsamps-1]'/samprate;
dt = 1/samprate;
T = tt(end) + dt; % total time length

if resamprate<samprate % else, will be upsampling - just interp;
    %taper ends
    for ii = 1:size(dat,2)
        dat_win(:,ii) = flat_hanning_win(tt,dat(:,ii),0,tt(end),T/200);
    end
    % apply low-pass filter 
    dat_win_filt = filt_quick( dat_win,0,resamprate/2,1./samprate);
else
    dat_win_filt = dat;
end

% down-sampled time axis
dt_down = 1/resamprate;
tt_down = [0:dt_down:(T-dt_down)]';
% interpolate to downsample
dat_down = interp1(tt,dat_win_filt,tt_down,'linear','extrap');
% here we extrapolate linearly outside of tt - in the case that we are
% upsampling. 

if sum(tt_down>max(tt))*dt_down > dt % shouldn't be more to extrapolate than a single original dt
    error('seem to be extrpolating too far')
end

% plots to check
% plot(tt,dat,'k', tt, dat_win,'r',tt,dat_win_filt,'g')
% plot(tt,dat,'k',tt_down,dat_down,'r')

end

