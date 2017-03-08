function [ visc_melt,visc_hk ] = melt_visc_reduce_TH09( phi, visc_nomelt )
% [ visc_melt,visc_hk ] = melt_visc_reduce_TH09( phi, visc_nomelt )
% 
%   Function to apply the contiguity model of Takei and Holtzman 2009a
% 
%  phi is the porosity, i.e volumetric melt fraction 
% 
% eta_wmelt = 0.2 * contig^2 * eta_nomelt
% contig = 1 - A * sqrt(phi)


if nargin < 2 || isempty(visc_nomelt)
    visc_nomelt = 1;
end

A = 2.3;

contig = 1 - A * sqrt(phi);

visc_melt = 0.2 * contig.^2 * visc_nomelt;

lambda = 26; % for diffusion creep from Mei et al. 2002
visc_hk = 0.2 * exp(-lambda*phi) * visc_nomelt;
% fprintf('For reference, the hirth/kohlstedt value is %.2f greater\n',visc_hk/visc_melt)

visc_melt(phi==0) = visc_nomelt;
visc_hk(phi==0) = visc_nomelt;


end

