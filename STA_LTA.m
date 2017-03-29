function [ t0,stalta ] = STA_LTA( dat,dt,tt )
% [ t0,ltasta ] = STA_LTA( dat,dt,[threshold] )

if nargin < 3 || isempty(tt)
    tt = dt*[0:size(dat,1)-1]';
end


lt_w = 5; % "long-term" running av. window 
st_w = 1; % "short-term" running av. window 

% pad data so that END of lt_w and st_w coincide
% square data to get power 
lt_dat = [zeros(0.5*(lt_w-st_w)./dt,1);dat.^2];
st_dat = [dat.^2;zeros(0.5*(lt_w-st_w)./dt,1)];

lta = moving_average(lt_dat,lt_w/dt);
sta = moving_average(st_dat,st_w/dt);
lta(lta<0.01*rms(lta)) = nan;

stalta = sta./lta;

% figure(1);clf 
% subplot(311);plot(dat);
% subplot(312);hold on;plot(lta);plot(sta);
% subplot(313);plot(sta./lta);

[~,ind] = maxab(stalta);
t0 = tt(ind);

end

