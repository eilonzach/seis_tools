function [ ott ] = migrate_PorS_conv( itt,model,rayp,rayp_ref,ph )
%[ ott ] = migrate_PorS_conv( itt,model,rayp,rayp_ref,ph )
% 
% Function to migrate a time series by scaling the time axes to the time
% axis that would be recovered for a ray with reference ray parameter. If
% the reference rayp is zero, then the output is a time series of depths. 
% 
% INPUTS:
%  itt   - vector of times corresponding to the data
%  model - velocity model; structure with fields 'z','VS', and 'VP'
%  rayp  - ray parameter of the data, in units of s/km
%  rayp_ref - reference ray parameter (s/km). If this is ==0, then the
%             output is a series of depths (in km), not times
%  ph    - phase and conversion, 'Ps' or 'Sp' 
% 

zz = model.z(:);
vs = model.VS(:);
vp = model.VP(:);

%%% Can work this code into the function, as an option, for doing the Earth
%%% flattening transformation to handle spherical Earth. Makes a small
%%% change to results. 
% a = 6371; 
% r = a - zz; 
% zzf = -a .* log(r./a); 
% vsf = a./r .* vs; 
% vpf = a./r .* vp; 
% 
% zz = zzf; 
% vs = vsf; 
% vp = vpf; 
%%%

if strcmp(ph,'Ps')
    tt_mig = diff(zz).*( sqrt(midpts(1./vs./vs) - rayp^2) - sqrt(midpts(1./vp./vp) - rayp^2) );
    if rayp_ref~=0
        tt_mig_ref = diff(zz).*( sqrt(midpts(1./vs./vs) - rayp_ref^2) - sqrt(midpts(1./vp./vp) - rayp_ref^2) );
    end
elseif strcmp(ph,'Sp')
    tt_mig = diff(zz).*( sqrt(midpts(1./vp./vp) - rayp^2) - sqrt(midpts(1./vs./vs) - rayp^2) );
    if rayp_ref~=0
        tt_mig_ref = diff(zz).*( sqrt(midpts(1./vp./vp) - rayp_ref^2) - sqrt(midpts(1./vs./vs) - rayp_ref^2) );
    end
end

% tt_mig = cumtrapz([0;tt_mig]); % Trapz gives a max of about 0.1 difference over 40 s versus cumsum. Not necessary. Comment is left in case further exploration is desired. 
tt_mig = cumsum([0;tt_mig]); 
if rayp_ref~=0  
%     tt_mig_ref = cumtrapz([0;tt_mig_ref]);
    tt_mig_ref = cumsum([0;tt_mig_ref]);
end

ind = [1:length(tt_mig)]; 

% cut out repetitions in the vectors
reps = tt_mig(1:end-1) == tt_mig(2:end); 
reps = [reps; false]; % Make sure reps has same size as zz, vs, and vp. 

% cut out nans (i.e. inhomogeneous) from both vectors
if rayp_ref~=0
    inhom = (imag(tt_mig)~=0) | (imag(tt_mig_ref)~=0);
else
    inhom = (imag(tt_mig)~=0);
end

ind( or(inhom, reps) ) = [];

if rayp_ref~=0
    ott = interp1(tt_mig(ind),tt_mig_ref(ind),itt,'linear',nan);
else
    ott = interp1(tt_mig(ind),zz(ind),itt,'linear',nan);    
end
    
end

