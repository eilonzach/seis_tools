function [K,G] = KG_from_VpVsrho(vp,vs,rho)
% [K,G] = KG_from_VpVsrho(vp,vs,rho)
%
% Function to compute shear and bulk moduli from velocities
% Units of moduli in GPa
% Units of velocity in km/s, density in kg/m^3. 
% Will attempt to recognise and adjust if velocities in m/s or rho in g/cc

if vp > 1e2
    fprintf('Looks like Vp is in m/s, adjusting units')
    vp = vp/1e3;
end
if vs > 1e2
    fprintf('Looks like Vs is in m/s, adjusting units')
    vs = vs/1e3;
end
if rho < 1e2
    fprintf('Looks like rho is in g/cc, adjusting units')
    rho = rho*1e3;
end

% put in SI units
vp = vp*1e3;
vs = vs*1e3;

G = vs.^2 * rho;
K = vp.^2 * rho - (4/3)*G;
G = G/1e9;
K = K/1e9;


end