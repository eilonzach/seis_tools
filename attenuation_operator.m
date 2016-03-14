function [ Dwt ] = attenuation_operator( Q0,c0,L,w,alpha,opt )
%[ Dwt ] = attenuation_operator( Q,L,omega,alpha )
%   Function to calculate the attenuation operator, the complex function in
%   frequency space that applies a given attenuation - i.e. multiply this
%   by the fourier transform of a signal and ifft to get the attenuated
%   signal.
% 
% INPUTS:
%  Q0      - the Quality factor (in the case of f-dependent Q, the Q0 at 1Hz)
%  L       - the distance travelled, in m
%  c0      - the reference phase velocity (at 1Hz), in m/s
%  w       - the vector of angular frequencies - this can be just the
%             positive, real frequencies or the freqs associated with fft
%  alpha   - [optional, default=0] the frequency dependence of q
%  opt     - [optional, default='zph'] choice to have zero-phase 'zph' or
%             minimum phase 'mph'

% OUTPUT:
%   Dwt - the attenuation operator, a complex function of length same as
%         omega

if nargin<5
    alpha = 0;
end
if nargin<6
    opt = 'zph';
end

qinv = 1/Q0;

if alpha==0
    qinv_w = qinv;
    c_w = c0 * (1 + (1/(Q0*pi))*log(abs(w)./2./pi)); 
    % c_max = c0 * (1 + (1/(Q*pi))*log(max(w)./2./pi));
    c_inf = c0 * (1 + (1/Q0)); % artificial - actually goes non causal
elseif alpha~=0
    qinv_w = qinv.*abs(w/2/pi).^(-alpha);
    qinv_w(1) = 2.5*qinv_w(2)-qinv_w(3); % fudge to account for zero-freq.
    % USING Karato 1993 eq(5) or Minster and Anderson 1981 eq(43) where c0
    % is defined as the anharmonic velocity (i.e. w_ref = Inf)
    c_w = c0 * (1 - 0.5*cot(0.5*pi*alpha)*qinv_w); 
    c_inf = c0; % precise, because the f-dependent term is causal
    % USING Lekic 2009 eq(11) with reference freq w=1;
%     c_w = c0 * (1 + (0.5*cot(0.5*pi*alpha)*qinv)*(1 - abs(w/2/pi).^-alpha) );
% 	c_inf = c0 * (1 + (0.5*cot(0.5*pi*alpha)*qinv) );
end




c_w = c_w(:);


if strcmp(opt,'zph')
%% zero-phase attenuation operators (Fang & Muller 1991)
% equation 1 - Delayed attenuator operator in the f domain
Dwt = exp( -0.5*abs(w)*L.*qinv_w./c0 ) .* exp( -1i*w.*L.*(1./c_w - 1./c_inf) );
elseif strcmp(opt,'mph')
% equation 5 - complete minimum-phase attenuation operator
Dwt = exp( -0.5*abs(w)*(L/c0) .* ( qinv_w + 1i*imag(hilbert(qinv_w)) ) );
end



end

