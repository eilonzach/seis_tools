function [ t0,stalta ] = STA_LTA( dat,dt,tt,threshold,lt_w,st_w )
% [ t0,ltasta ] = STA_LTA( dat,dt,tt,threshold,lt_w,st_w )

if nargin < 3 || isempty(tt)
    tt = dt*[0:size(dat,1)-1]';
end
if nargin < 4 || isempty(threshold)
    threshold=1.8;
end
if nargin < 5 || isempty(lt_w)
    lt_w = 15; % "long-term" running av. window 
end
if nargin < 6 || isempty(st_w)
    st_w = 2; % "short-term" running av. window 
end


% square data to get power 
lt_dat = dat.^2;
lta = moving_average(lt_dat,lt_w/dt);

% subtract a moving average from the data in the sttime series
st_dat = [dat - moving_average(dat,lt_w/dt)].^2; 
st_dat = [dat].^2; 
sta = moving_average(st_dat,st_w/dt);

lta(lta<0.01*rms(lta)) = nan;
sta(sta<0.1*max(sta)) = nan;

stalta = sta./lta;
% stalta_rms = sqrt(nanmean(stalta.^2))
% stalta(stalta<nanrms(stalta)) = nan;

% figure(1);clf, set(gcf,'pos',[186 396 560 645])
% subplot(411);plot(tt,dat),xlim([tt(1),tt(end)]), 
% subplot(412);hold on; plot(lt_dat); plot(st_dat);xlim([0,length(stalta)])
% subplot(413);hold on;plot(lta);plot(sta);xlim([0,length(stalta)])
% subplot(414);plot(stalta);xlim([0,length(stalta)])

ind = find(stalta>threshold,1,'first');
t0 = tt(ind) ;
if isempty(t0), t0=nan; end

end

