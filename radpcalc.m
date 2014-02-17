function [ R ] = radpcalc(str,dip,rak,faz,inc)
%[ R ] = radpcalc(str,dip,rak,faz,inc)
%  For an earthquake with geometry defined by strike, dip and rake, this
%  script takes input vectors of azimuth and incidence angle and computes
%  the radiation amplitudes for each ray leaving the source
%
% INPUTS
% str   = fault strike, degrees from north
% dip   = fault dip, degrees from horizontal
% rak   = fault rake, +90 is pure thrust, 0 is left lateral for a ss
% faz   = vector of leaving azimuths of rays, degrees from north
% inc   = vector of leaving angles of rays, degrees from vertical
%   n.b. faz and inc can be vectors (length N) or scalars, if both are
%   vectors, they have to be the same length
% 
% OUTPUT
%   R   = 3xN matrix, each row is relative (scaled to 1) amplitude of 
%       [ P SV SH ]

if isrow(faz); faz=faz'; end
if isrow(inc); inc=inc'; end

str = d2r(str)*ones(size(faz));
dip = d2r(dip);
rak = d2r(rak);
faz = d2r(faz); 
i = d2r(inc);

% Equations from Stein and  Wysession p233 
sR = sin(rak)*sin(dip)*cos(dip);
qR = sin(rak)*cos(2*dip)*sin(str-faz) + cos(rak)*cos(dip)*cos(str-faz);
pR = cos(rak)*sin(dip)*sin(2*(str-faz)) - sin(rak)*sin(dip)*cos(dip)*cos(2*(str-faz));
qL = -cos(rak)*cos(dip)*sin(str-faz) + sin(rak)*cos(2*dip)*cos(str-faz);
pL = sin(rak)*sin(dip)*cos(dip)*sin(2*(str-faz)) + cos(rak)*sin(dip)*cos(2*(str-faz));

Rp = sR.*(3*(cos(i).^2) - 1) - qR.*sin(2*i) - pR.*(sin(i).^2);
Rsv = 1.5*sR.*sin(2*i) + qR.*cos(2*i) + 0.5*pR.*sin(2*i);
Rsh = -qL.*cos(i) - pL.*sin(i);

R = [Rp, Rsv, Rsh];


end

