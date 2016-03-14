function [ Mw_or_Mo ] = Mw_calc( option, Mo_or_Mw, unit )
% [ Mw_or_Mo ] = Mw_calc( option, Mo_or_Mw, unit )
%   function to calculate moment magnitude from mantissa and exponent of
%   moment, or vice versa!
%   can do units of dynecm  (default) or Nm
% USAGE:
% Mw = Mw_calc('Mo2Mw',1.0e25)
% Mw = Mw_calc('Mo2Mw',1.0e25,unit) - where unit = 'dynecm'/'Nm'
% Mo = Mw_calc('Mw2Mo',6.2)
% Mo = Mw_calc('Mw2Mo',6.2, unit)
% 
% Z. Eilon 2015

if nargin < 3,
    unit = 'dynecm';
end

if      strcmp(unit,'dynecm'), f = 1;
elseif  strcmp(unit,'Nm'),     f = 1e7;
end

if strcmp(option,'Mo2Mw') | strcmp(option,'M02Mw')
    Mo = Mo_or_Mw;
    Mw = (2/3)*log10(Mo*f) - 10.7;
    Mw_or_Mo = Mw;
elseif strcmp(option,'Mw2Mo') | strcmp(option,'Mw2M0')
    Mw = Mo_or_Mw;
    Mo = 10^(1.5*(Mw + 10.7))/f;
    Mw_or_Mo = Mo;
else
    error('Option not recognised')
end
    


end

