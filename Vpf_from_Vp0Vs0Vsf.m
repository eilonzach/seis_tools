function [ Vp_f ] = Vpf_from_Vp0Vs0Vsf( Vp_0,Vs_0,Vs_f )
%[ Vp_f ] = Vpf_from_Vp0Vs0Vsf( Vp_0,Vs_0,Vs_f )
%   Assuming no attenuation in bulk (i.e. Q_kappa infinite) calculate the
%   change in Vp due to relaxation of the shear modulus, computed using the
%   unrelaxed Vp, and the anelastic + unrelaxed Vs


Vp_f = sqrt( Vp_0.^2 + (4/3)*(Vs_f.^2 - Vs_0.^2) );



end

