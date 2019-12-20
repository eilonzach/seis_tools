function [ dat_t,tt ] = migrate_depth2time(dat_z,model,rayp,ph,dt )
% [ dat_t,tt ] = migrate_depth2time(dat_z,model,rayp,ph,dt )
% 
% Function to migrate a time series by scaling the time axes to the time
% axis that would be recovered for a ray with reference ray parameter. If
% the reference rayp is zero, then the output is a time series of depths. 
% 
% INPUTS:
%  dat_z - vector of data as a function of depth - same depths as in model
%  model - velocity model; structure with fields 'z','VS', and 'VP'
%  rayp  - ray parameter of the desired output data, in units of s/km
%  ph    - phase and conversion, 'Ps' or 'Sp' 
%  dt    - dt for output time series (1/samprate) 
% 

if rayp > 1
    fprintf('\n !! rayp seems to be in s/deg in migrate_depth2time...\n       ...automatically converting to s/km\n\n');
    rayp = rayp_sdeg2skm(rayp);
end

zz = model.z(:);
vs = model.VS(:);
vp = model.VP(:);

tt_mig = diff(zz).*( sqrt(midpts(1./vs./vs) - rayp^2) - sqrt(midpts(1./vp./vp) - rayp^2) );
tt_mig = cumsum([0;tt_mig]);


% regularly spaced time vector
tt = [0:dt:max(tt_mig)]';

% interp to get time-sampled data vector
dat_t = interp1(tt_mig,dat_z,tt,'linear',nan);

if strcmp(ph,'Ps')

elseif strcmp(ph,'Sp')
    tt = -tt; % flip time axis - deeper structure is earlier and earlier...
end


% % cut out repetitions in the vectors
% reps = find(diff(tt_mig)==0);
% ind(reps)= []; %#ok<FNDSB>
% % cut out nans (i.e. inhomogeneous) from both vectors
% if rayp_ref~=0
%     inhom = (imag(tt_mig)~=0) | (imag(tt_mig_ref)~=0);
% else
%     inhom = (imag(tt_mig)~=0);
% end
% ind(inhom) = [];


end

