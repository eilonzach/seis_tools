function [ R ] = radpcalc(str,dip,rak,esaz,inc)
%[ R ] = radpcalc(str,dip,rak,faz,inc)
%  For an earthquake with geometry defined by strike, dip and rake, this
%  script takes input vectors of azimuth and incidence angle and computes
%  the radiation amplitudes for each ray leaving the source
%
% INPUTS
% str   = fault strike, degrees from north
% dip   = fault dip, degrees from horizontal
% rak   = fault rake, +90 is pure thrust, 0 is left lateral for a ss
% esaz   = vector of leaving azimuths of rays, degrees from north
% inc   = vector of leaving angles of rays, degrees from vertical
%   n.b. faz and inc can be vectors (length N) or scalars, if both are
%   vectors, they have to be the same length
% 
% OUTPUT
%   R   = 3xN matrix, each row is relative (scaled to 1) amplitude of 
%       [ P SV SH ]

esaz=esaz(:);
inc=inc(:);

str = d2r(str)*ones(size(esaz));
dip = d2r(dip);
rak = d2r(rak);
esaz = d2r(esaz); 
i = d2r(inc);

% % Equations from Stein and  Wysession p233 
% sR = sin(rak)*sin(dip)*cos(dip);
% qR = sin(rak)*cos(2*dip)*sin(str-esaz) + cos(rak)*cos(dip)*cos(str-esaz);
% pR = cos(rak)*sin(dip)*sin(2*(str-esaz)) - sin(rak)*sin(dip)*cos(dip)*cos(2*(str-esaz));
% qL = -cos(rak)*cos(dip)*sin(str-esaz) + sin(rak)*cos(2*dip)*cos(str-esaz);
% pL = sin(rak)*sin(dip)*cos(dip)*sin(2*(str-esaz)) + cos(rak)*sin(dip)*cos(2*(str-esaz));
% 
% Rp = sR.*(3*(cos(i).^2) - 1) - qR.*sin(2*i) - pR.*(sin(i).^2);
% Rsv = 1.5*sR.*sin(2*i) + qR.*cos(2*i) + 0.5*pR.*sin(2*i);
% Rsh = -qL.*cos(i) - pL.*sin(i);

% Equations from Aki and Richards p108-109
dphi = esaz-str;
cR = cos(rak);
sR = sin(rak);
cD = cos(dip);
sD = sin(dip);
cI = cos(i);
sI = sin(i);

s2D = sin(2*dip);
c2D = cos(2*dip);
s2I = sin(2*i);
c2I = cos(2*i);

Rp =    cR.*sD.*sI.*sI.*sin(2*dphi)                     ...
      - cR.*cD.*s2I.*cos(dphi)                           ...
      + sR.*s2D.*(cI.*cI - sI.*sI.*sin(dphi).*sin(dphi))...
      + sR.*c2D.*s2I.*sin(dphi);
  
Rsv =   sR.*c2D.*c2I.*sin(dphi)                         ...
      - cR.*cD.*c2I.*cos(dphi)                          ...
      + 0.5*cR.*sD.*s2I.*sin(2*dphi)                    ...
      - 0.5*sR.*s2D.*s2I.*(1 + sin(dphi).*sin(dphi));
  
Rsh =   cR.*cD.*cI.*sin(dphi)                           ...
      + cR.*sD.*sI.*cos(2*dphi)                         ...
      + sR.*c2D.*cI.*cos(dphi)                          ...
      - 0.5*sR.*s2D.*sI.*sin(2*dphi);




R = [Rp, Rsv, Rsh];


end

