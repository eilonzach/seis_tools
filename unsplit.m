function [ Ru,Tu ] = unsplit( Rs, Ts, phi, dT, samprate, inipol )
%[ Ru,Tu ] = unsplit( Rs, Ts, phi, dT, samprate, inipol )
%   function to perform the reverse splitting operation, attempting to undo
%   the effect of simple, one-stage splitting. If phi and dT are correct,
%   this should re-linearise elliptical particle motion and result in zero
%   tangential (assuming SKS) energy.
%
%  INPUTS
%   Rs      - Nx1 vector of split radial component
%   Ts      - Nx1 vector of split tangential component
%   phi     - azimuth of fast axis of anis (degrees)
%   dT      - splitting time (s)
%   samprate- sample rate of signal
%   inipol  - initial polarisation of Rs - i.e. forward azimuth of wave
%
%  OUTPUTS
%   Rs      - Nx1 vector of unsplit radial component
%   Ts      - Nx1 vector of unsplit tangential component

RR = Rs;
TT = Ts;

for ii = 1:length(phi)
daz = phi(ii)-inipol;

[ datFS ] = rot_traces( [RR,TT], daz);
sampshift = round(dT(ii).*samprate./2);
datF = [zeros(sampshift,1);datFS(1:end-sampshift,1)];
datS = [datFS(sampshift+1:end,2);zeros(sampshift,1)];
[ datRT2 ] = rot_traces( [datF,datS], -daz);
RR = datRT2(:,1);
TT = datRT2(:,2);
end

Ru = RR;
Tu = TT;


%% Junk used while testing...
% samprate= 50;
% phi     = 143;
% dT      = 0.6;
% 
% %  ===== for testing ====
% 
% samplength=100;
% impulselength=20;
% amp=10;
% startfraction=0.4;
% waveform='gauss';
% dt = 1./samprate;
% inipol = 100;
% 
% xx = synthtrace(samplength,impulselength,amp,dt,waveform,startfraction);
% tt = [0:dt:samplength-dt];
% [odat1,odat2] = split1layer(xx,zeros(size(xx)),50,inipol,phi,dT,'RT');
% 
% 
% Rs = odat1;
% Ts = odat2;
% 


% %  ======= post hoc testing ========
% 
% figure(1), clf, 
% subplot(3,1,1)
% plot(tt,xx)
% subplot(3,1,2)
% plot(tt,[odat1,odat2])
% subplot(3,1,3)
% plot(tt,datRT2)

end

