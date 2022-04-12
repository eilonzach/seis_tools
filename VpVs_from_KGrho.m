function [vp,vs] = VpVs_from_KGrho(K,G,rho)
% [vp,vs] = VpVs_from_KGrho(K,G,rho);
%
% Function to compute velocities from elastic moduli and density
% Units of velocity in km/s
% Units of moduli in GPa, density in kg/m^3. Will attempt to recognise and
% adjust if values in Pa or g/cc.

if K > 1e3
    fprintf('Looks like K is in Pa, adjusting units')
    K = K/1e9;
end
if G > 1e3
    fprintf('Looks like G is in Pa, adjusting units')
    G = G/1e9;
end
if rho < 1e2
    fprintf('Looks like rho is in g/cc, adjusting units')
    rho = rho*1e3;
end

% put in SI units
K = K*1e9;
G = G*1e9;

vp = sqrt( (K + (4/3)*G)/rho ); vp = vp/1e3;
vs = sqrt( G/rho );             vs = vs/1e3;


end