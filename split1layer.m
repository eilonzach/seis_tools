function [odat1,odat2] = split1layer(idat1,idat2,samprate,inipol,phi,dtime,coords)
%function [odat1,odat2] = split1layer(idat1,idat2,samprate,inipol,phi,dtime,coords)
% 
% Function to take an input waveform with a given polarisation and split it
% into fast and slow waves, before recombining onto two orthogonal channels
% (N&E or R&T)
%
% INPUTS:
% dat1     - time series of data on channel 1
% dat2     - time series of data on channel 2, orie 90? CLOCKWISE of ch1
%             can just be zeros, if want to assume that initial data is
%             linearly polarised
%            N.B. ensure there is more than dtime worth of padding on either
%             end of the time series, or useful signal might get cut off
%             Ideally, they would be pre-windowed, with zeros added
%
% samprate - sample rate of input time series (samples per second)
% inipol   - initial polarisation of ch1 (degrees clockwise from North)
% phi      - azimuth of fast axis of anisotropy (degrees clockwise from N)
% dtime    - time delay between fast and slow exiting waves 
%             each will be offset from initial arrival by ±dtime/2
%          - as a consequence, should be divisible by dt
% coords   - option for output: 
%               'NE' will make odat1 North and odat2 East
%               'RT' will make odat1 parallel to inipol and odat2 90?
%               clockwise from it
% all angles in DEGREES
% 
% OUTPUTS:
% odat1     - output time series along chan 1
% odat2     - output time series along chan 2, orie 90? CLOCKWISE of ch1

if nargin < 7 % no 'coords' option specified
    coords = 'NE'; % assume NE coord system
end

%% Rotate initial data into fast-slow coord system
daz=phi-inipol; % clockwise angle from inipol to phi
caz=cosd(daz);
saz=sind(daz);
rmat=[caz saz; -saz caz];
% make sure time series are in column vectors
if size(idat1,1)==1; idat1=idat1'; end
if size(idat2,1)==1; idat2=idat2'; end
% do rotation
ff =  caz*idat1 +  saz*idat2; % fast axis
ss = -saz*idat1 +  caz*idat2; % slow axis +ive to right, looking in +ive fast direction
% ff = (rmat(1,:)*[idat1';idat2'])'; % fast azis
% ss = (rmat(2,:)*[idat1';idat2'])'; % slow axis +ive to right, looking in +ive fast direction


%% apply time shifts to both traces
% have to add buffer - added sinusoidal decay to zero at both ends
sampshift = round(samprate*dtime/2); % number of samples to be removed and added
% nsamps=length(idat1);

ff = [ff(sampshift+1:end);ff(end)*cos(pi*[0:sampshift-1]'/(2*sampshift))];
ss = [ss(1)*sin(pi*[0:sampshift-1]'/(2*sampshift));ss(1:end-sampshift)];

if strcmp(coords,'NE')==1
    daz=-phi; % rotate all back to NE
elseif strcmp(coords,'RT')==1
    daz=-daz; % rotate back to original polarisation
end
caz=cosd(daz);
saz=sind(daz);
odat1 =  caz*ff +  saz*ss; 
odat2 = -saz*ff +  caz*ss; %orie 90? CLOCKWISE of odat1

end




