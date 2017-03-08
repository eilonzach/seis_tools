
function test_ETMTM()
    test_ETMTM_nick()
end

function test_ETMTM_nick()
    %A test routine for ETMTM
    %
    clf;
    for lowT = [2.0 4.0 8.0];

        npts=2000;
        
        %Initial traces are white noise
        P=rand(npts,1)'.*0.01;
        D=rand(npts,1)'.*0.01;
        
        %Add three spikes of varying amplitude
        P(1750)= 1;
        D(1750)= 0.5;
        D(1500)= -0.4;
        D(1400)= 0.3;

        dt = 0.05;

        %lowT=2.0;
        highT=100.0;

        lf=1.0/highT;
        hf=1.0/lowT;

        P=bpfilt(P,dt,lf,hf);
        D=bpfilt(D,dt,lf,hf);

        %Apply a taper
        T_Wt            = 20;               % Taper width (sec)
        T_W         	= ceil(T_Wt/dt);
        T             	= ones(1,length(P));
        T(1:T_W+1)    	= (1-cos((0:T_W)*pi/T_W))/2;
        T(end-T_W:length(P))	= (1+cos((0:T_W)*pi/T_W))/2;
        
        P_T         = P.*T;   	% Parent tapered
        D_T         = D.*T;    	% Daughter tapered

        TB = 4;
        NT = 7;

        %RF_Time_Tap = ETMTM(P_T,D_T,TB,NT,t0,t1,'Sp','synth',dt);
        [Time, RF_Time] = ETMTM(P_T,D_T,TB,NT,'synth',dt);

        plot(Time,RF_Time)
        xlabel('Time (s)')
        ylabel('RF Amplitude')
        hold on
    
    end
end

function [y] = bpfilt(x,dt,lf,hf)

% ********* Function Description *********
%
% Bandpass filter a time seriers.
%
% [Y] = bpfilt(X,DT,LF,HF)
%
% Take a time series, X, sampled at DT and
% filter it with a 2nd order, 2 pass
% butterworth filter between frequencies
% LF and HF. If X is a matrix, this will
% filter the individual rows of X.
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Taken from code written by Michael  *
% *  Bostock and Stephane Rondenay, and  *
% *  used by Kate Rychert.               *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

nyq     = 0.5/dt;           % Nyquist Frequency
wn      = [lf/nyq,hf/nyq];
[b,a]   = butter(2,wn);
for ix=1:length(x(:,1))
    y(ix,:) = filtfilt(b,a,double(x(ix,:)));  % Edited by Ved because waveforms
    % are single not double, but filtfilt
    % requires double
end

end
