function [ Vsv,Vsh ] = VsvVsh_from_VsXi( Vs,xi )
%[ Vsv,Vsh ] = VsvVsh_from_VsXi( Vs,Xi )
%   Function to calculate Vsv and Vsh from the voigt average velocity (Vs)
%   and the value of xi, which describes the radial anisotropy, where
% 
%   Vs^2 = (Vsh^2 + 2Vsv^2)/3
%   xi = Vsh^2/Vsv^2 (=N/L)


Vsv = Vs .* sqrt(3./(xi+2));
Vsh = Vs .* sqrt(3.*xi./(xi + 2));


end

