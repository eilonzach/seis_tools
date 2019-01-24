function [ Vs_vt, Vp_vt ] = voigtav( Vsh,Vsv,Vph,Vpv,eta )
%[ Vs_av, Vp_av ] = voigtav( Vsh,Vsv,Vph,Vpv,eta )
%   Calculate Voigt average velocities from vertically and horizontally
%   polarised wavespeeds, plus eta...
% 
% NB removed rho from these equations - cancels...
% 
% If using approx expressions (eta=1)
%   Vs_vt = sqrt( (  Vsh.^2 + 2*Vsv.^2)/3 );
%   Vp_vt = sqrt( (4*Vph.^2 +   Vpv.^2)/5 );


if nargin > 2

c11 = Vph.^2;
c22 = c11;
c33 = Vpv.^2;
c44 = Vsv.^2;
c55 = c44;
c66 = Vsh.^2;
c12 = c11 - 2*c66;
c13 = eta.*(c11 - 2*c44);
c23 = c13;

K = (c11 + c22 + c33 + 2*(c12 + c13 + c23))/9;
G = (c11 + c22 + c33 + 3*(c44 + c55 + c66) - c12 - c13 - c23)/15;

Vs_vt = real(sqrt(G));
Vp_vt = real(sqrt(K + (4/3)*G));

else % just Vsh and Vsv ==> use approx expressions (i.e. assume eta=1)
    
Vs_vt = real(sqrt( (  Vsh.^2 + 2*Vsv.^2)/3 ));
Vp_vt = real(sqrt( (4*Vph.^2 +   Vpv.^2)/5 ));

end

end

