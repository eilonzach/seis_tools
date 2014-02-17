function [odat1,odat2] = splitNlayer(idat1,idat2,samprate,inipol,phis,dtimes,coords)
% function [odat1,odat2] = splitNlayer(idat1,idat2,samprate,inipol,phis,dtimes,coords)
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
% inipol   - initial polarisation of ch1 (degrees clockwise from North)%
% coords   - option for output: 
%               'NE' will make odat1 North and odat2 East
%               'RT' will make odat1 parallel to inipol and odat2 90?
%               clockwise from it
% all angles in DEGREES
% 
% LAYERs (first in the vectors are deeper, i.e. first encountered)
% phis      - vector with azimuths of fast axis of anisotropy (degrees 
%             clockwise from N)
% dtimes    - vector with time delays between fast and slow exiting waves 
%             each will be offset from initial arrival by ±dtime/2
% 
% OUTPUTS:
% odat1     - output time series along chan 1
% odat2     - output time series along chan 2, orie 90? CLOCKWISE of ch1

if nargin < 7 % no 'coords' option specified
    coords = 'NE'; % assume NE coord system
end

N = length(phis);

for ii=1:N-1
    [rdat,tdat] = split1layer(idat1,idat2,samprate,inipol,phis(ii),dtimes(ii),'RT');
    idat1=rdat;
    idat2=tdat;
end

[odat1,odat2] = split1layer(idat1,idat2,samprate,inipol,phis(end),dtimes(end),coords);



% end

