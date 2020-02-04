function [ Vpv,Vph ] = VpvVph_from_VpPhi( Vp,phi )
% [ Vpv,Vph ] = VpvVph_from_VpPhi( Vp,phi ) 
%   Function to calculate Vpv and Vph from the voigt average velocity (Vp)
%   and the value of phi, which describes the radial anisotropy, where
% 
%   Vp^2 = (4Vph^2 + Vpv^2)/5
%   phi = Vpv^2/Vph2 (=C/A)


Vph = Vp .* sqrt(5./(phi+4));
Vpv = Vp .* sqrt(5.*phi./(phi + 4));


end

