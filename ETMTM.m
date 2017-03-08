
function[Time, RF_Time] = ETMTM(P,D,Ps_or_Sp,TB,NT,tag,dt,win_len,Poverlap)
    % ETMTM(P,D,Ps_or_Sp,TB,NT,tag,dt,win_len,Poverlap)
    %
    %ETMTM routine originally coded up by Ved Lekic,
    %  updated by Nick Mancinelli in 2016.
    %  slightly edited by Z. Eilon (Nov 2016)
    %   
    %
    %
    % Inputs:
    %   P:  Parent waveform - the parent phase should peak at about
    %           (win_len / 2) s before the end (for S-to-P)
    %   D:  Daughter waveform 
    %   Ps_or_Sp: option for 'Ps' or 'Sp'
    %   TB: Period * Bandwith (e.g., 3)
    %   NT: Number of tapers (e.g., 7)
    %   tag: data or synth
    %   dt: delta t
    %   win_len: length of each window in seconds (e.g., 100)
    %   Poverlap: fractional window overlap (e.g., 0.90)
    %
    % Returns:
    %   RF_Time: a receiver function in time
    %
    % Findings: when P arrival is not in the center of the window, the
    % amplitudes are not unity at the beginning and decreasing from there on.
    % Instead they peak at the time shift which corresponds to the middle index
    % in the P time window.

    % As your TB
    % increases, the frequency smearing gets worse, which means that the RFs
    % degrate at shorter and shorter lag times. Therefore, as you increase TB,
    % you should also increase Poverlap.

    %TB = 4; NT = 7; %choise of TB = 4, NT = 3 is supposed to be optimal
    %t0 = -5; t1 = max(time);
    %        function [RF_Time] = MTMDecon_VedOptimized(P,D,TB,NT,t0,t1,Faza)
    % Ved wrote MTM for MATLAB, which has the added advantage of
    % finding the optimal damping parameter.
    % TB  = time bandwidth product (usually between 2 and 4)
    % NT  = number of tapers to use, has to be <= 2*TB-1

    t1              = dt*(length(P)-1);
    t0              = -5;
    
    % For Sp we have to flip time axis -- not for Ps
    if strcmp(Ps_or_Sp,'Sp'), 
        D = fliplr(D); P = fliplr(P);
    elseif strcmp(Ps_or_Sp,'Ps')
        %:)
    end

    % Length of moving time window in seconds
    % win_len = 50;
    Nwin = round(win_len/dt);

    % Fraction of overlap overlap between moving time windows. As your TB
    % increases, the frequency smearing gets worse, which means that the RFs
    % degrate at shorter and shorter lag times. Therefore, as you increase TB,
    % you should also increase Poverlap.
    %Poverlap = 0.90;
    nso=(1-Poverlap)*Nwin; nso=round(nso);

    npad=zeros(1,nso*1); D=[npad D npad]; P=[P npad npad];
    time = 0:dt:dt*(length(P)-1);

    % Create moving time windows and daughter/parent snippets
    starts = 1:round((1-Poverlap)*Nwin):length(P)-Nwin+1; nd=0;
    for j = 1:length(starts)
        tmp_times(j,1:Nwin) = time(starts(j):starts(j)+Nwin-1);
        if(j==1) % ASSUME THAT PARENT PHASE IS CENTERED IN FIRST WINDOW!
            Pwin = interp1(double(time),double(P),double(tmp_times(j,:)),'linear',0)';

        end
        Dwin(1:Nwin,j) = interp1(double(time),double(D),double(tmp_times(j,:)),'linear',0);

        ltp=win_len/5; % taper before deconvolving (important for synthetics)
        Dwin(1:Nwin,j)=taper(Dwin(1:Nwin,j)',tmp_times(j,1:Nwin),ltp,dt,...
            tmp_times(j,1+round((ltp+dt)/dt)),tmp_times(j,Nwin-round((ltp+dt)/dt)));
        nd1=length(find(Dwin(:,j)>1e-2)); if nd1>nd; nd=nd1; end
    end

    % Search range for optimal damping parameter alpha
    switch tag
        case 'data'
            alphas = logspace(-2,2,20)*var(D(round(end/4):3*round(end/4)))*length(P);
        case 'synth'
            alphas = logspace(-2,2,20)*var(D)*length(P);
    end
    % Figure out average times for each moving window
    t0_times = median(tmp_times,2);

    % Construct Slepians
    [E,~] = dpss(length(Pwin),TB);

    % Length of waveforms;
    nh = length(Pwin);
    %
    misfit = zeros(size(alphas)); magntd = zeros(size(alphas));
    % Now, calculate misfit and RF size for each alpha
    for kj = 1:length(alphas)

        % Now loop through different time windows for the daughter component
        for k = 1:size(tmp_times,1)
            % Create multitaper estimates
            for j = 1:NT
                tmp1 = fft(E(:,j).*Pwin,nh);
                tmp2 = fft(E(:,j).*Dwin(:,k),nh);
                if j==1
                    NUM = conj(tmp1).*tmp2;
                    DEN = conj(tmp1).*tmp1;
                else
                    NUM = NUM + conj(tmp1).*tmp2;
                    DEN = DEN + conj(tmp1).*tmp1;
                end
            end

            % Calculate optimal RF
            tmp = real(ifft(NUM./(DEN + alphas(kj))));

            % Filter and Normalize optimal RF
            nrm = max(real(ifft(DEN./(DEN + alphas(kj)))));

            tmp = reshape(fftshift(tmp./nrm),[1 nh]);

            % Time vector
            vrijeme = dt*(-0.5*(nh-1):1:0.5*(nh-1))+t0_times(k)-t0_times(1) ...
                -dt.*length(npad);

            RF_Time_win(:,k) = interp1(double(vrijeme),double(tmp),double(t0:dt:t1),'linear',NaN);
            %length_Dwin=time(Nwin)-time(1);
        end


        tmp = conv(nanmean(RF_Time_win,2),P); t0d=round(t0/dt);
        if size(tmp,1)<size(tmp,2); mfD=D; else mfD=D'; end
        misfit(kj) = nansum(abs(mfD - tmp(1-t0d:length(D)-t0d)));
        magntd(kj) = nansum(abs(tmp(1-t0d:length(D)-t0d)));
    end

    % Find optimal alpha
    [~,j2] = min((misfit./std(misfit)).^2+(magntd./std(magntd)).^2);

    % Now loop through different time windows for the daughter component
    for k = 1:size(tmp_times,1)
        % Create multitaper estimates
        for j = 1:NT
            tmp1 = fft(E(:,j).*Pwin,nh);
            tmp2 = fft(E(:,j).*Dwin(:,k),nh);
            if j==1
                NUM = conj(tmp1).*tmp2;
                DEN = conj(tmp1).*tmp1;
            else
                NUM = NUM + conj(tmp1).*tmp2;
                DEN = DEN + conj(tmp1).*tmp1;
            end
        end

        % Find optimal alpha
        % [junk,j2] = min(abs(misfit./std(misfit))+abs(magntd./std(magntd)));

        % Calculate optimal RF
        tmp = real(ifft(NUM./(DEN + alphas(j2))));

        % Filter and Normalize optimal RF
        nrm = max(real(ifft(DEN./(DEN + alphas(j2)))));

        tmp = reshape(fftshift(tmp./nrm),[1 nh]);

        % Time vector
        vrijeme = dt*(-0.5*(nh-1):1:0.5*(nh-1))+t0_times(k)-t0_times(1) ...
                -dt.*length(npad);

        RF_Time_win(:,k) = interp1(double(vrijeme),double(tmp),double(t0:dt:t1),'linear',0);
    end

    RF_Time = nanmean(RF_Time_win,2);

    %Flip if S to P
	if strcmp(Ps_or_Sp,'Sp'), 
           Time = -(t0:dt:t1);
    elseif strcmp(Ps_or_Sp,'Ps')
           Time = (t0:dt:t1);
    end

end


function [y] = taper(x,time,tpr,dt,t1,t2)
        
        % ********* Function Description *********
        %
        % TAPER  Taper a time series.
        %
        % TAPER(X,TIME,TPR,DT,T1,T2) takes time
        % series sampled at DT and tapers it with
        % a cosine taper TPR seconds long from
        % beginning point T1-TPR and with reverse
        % cosine taper from point T2 to point T2+
        % TPR. Points outside the range (T1-TPR,
        % T2+TPR) are zeroed. If T1/T2 is negative
        % then taper is not implemented at the
        % beginning/end. If X is an array of
        % seismograms, then the taper is applied
        % to each row of X.
        %
        %
        % ****************************************
        % *                                      *
        % *  Modified from Kate Rychert's        *
        % *  receiver function code - May 2008   *
        % *                                      *
        % *  Email: David_Abt@brown.edu          *
        % *                                      *
        % ****************************************
        
        nn      = length(x(1,:));
        nx      = length(x(:,1));
        taper   = ones(1,nn);
        it      = (0:fix(tpr/dt))*dt/tpr;
        ct      = 0.5-0.5*cos(pi*it);
        T1      = fix(time(1)/dt);             % Absolute sample point of first time step
        it1     = fix(t1/dt+1)-T1;
        it2     = fix(t2/dt+1)-T1;
        
        % Emily, 6th May 2013: problem with start time of phase being <100s, so
        % taper subscripts are negative.
        % Temporary workaround, set taper length to be shorter.
        % N.B. Only one event so far has had this issue!
        
        if t1>0
            if it1>fix(tpr/dt)
                taper(it1-fix(tpr/dt):it1)    = ct;
                taper(1:it1-fix(tpr/dt))    = zeros(size(1:it1-fix(tpr/dt)));
            else
                taper(1:it1) = ct(fix(tpr/dt)-it1+2:end);
                taper(1) = 0;
                disp('Bizarre taper!')
                
            end
        end
        if t2>0
            if t2>time(end)-tpr
                t2  = time(end)-tpr;
                it2    = fix(t2/dt)-T1;
            end
            taper(it2:it2+fix(tpr/dt))    = fliplr(ct);
            taper(it2+fix(tpr/dt):nn)    = zeros(size(taper(it2+fix(tpr/dt):nn)));
        end
        
        y = zeros(nx,nn);
        for ix=1:nx
            y(ix,:) = x(ix,:).*taper;
        end
        
end