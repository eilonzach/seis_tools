function [ Z, X, Y ] = radpplot(str,dip,rak,psvsh,fign,subp);
% [ Z, X, Y ] = radpplot(str,dip,rak,psvsh,fign,subp)
% This script plots the radiation pattern for the input fault geometry and
% the given phase. If no phase is given, does P
% 
%INPUTS
% str   = fault strike, degrees from north
% dip   = fault dip, degrees from horizontal
% rak   = fault rake, +90 is pure thrust, 0 is left lateral for a ss
% psvsh = p = 1,  sv = 2,  sh = 3 - if no input, does p
% fign  = figure number to do plots in 
% subp  = subplot indices - vector [y x i]
%
%OUTPUTS
% Z     = matrix of normalised amplitudes, (naz,nr,phase)
% X     = grid of X locations of Z points
% Y     = grid of Y locations of Z points
% 
% for plotting more on same axes, use conversion:
% r = r_max * (sin(inc/2)/sin(45))
% x = sin(az)*r';
% y = cos(az)*r';
% for this script, as is, rmax is 4
%
%   Z. Eilon    20 May 2013

phases = {'p','sv','sh'};
if nargin < 4 
    psvsh = 1;
end

if nargin < 5 
fign = 1;
end

if nargin < 6 
subp = [ 1 1 1 ];
end

% set up grid of radii and angles for the radiation pattern
Nr = 50;
Naz = 360;
rmax = 4;
r = linspace(0,rmax,Nr)'; % each column a different r
inc = 2*asin(sind(45)*r/rmax);
% for polar plot, r = r_max * (sin(inc/2)/sin(45))
%                   = r_max * (cos(dip/2) - sin(dip/2))
az = linspace(0,2*pi,Naz+1)'; % each row a different az
X = sin(az)*r';
Y = cos(az)*r';
Z = zeros(Naz+1,Nr,3); % each row is a different azimuth, column is different r, layers are P, SV,SH
for iz = 1:length(az)
R = radpcalc(str,dip,rak,r2d(az(iz)),r2d(inc));
Z(iz,:,:) = R;
end    
%% make nodal planes
% consider a set of 181 vectors spaced 1degree apart in the nodal plane,
% making an angle alpha with the direction of the strike, with positive
% alpha dipping to -Z. can work out the dip of these vectors 
%   D = asin(sin(alpha)*sin(dip1))
% and angle from the strike in the horiz plane:
%   phi = atan(tan(alpha)*cos(dip1))
% then x and y in the frame of the plane are just 
%   x_ = sin(phi)*r_
%   y_ = cos(phi)*r_
% where r_ = r_max * (cos(D/2) - sin(D/2))
% and a simple rotation by the strike is enough to get x_ and y_ to x and y
faz = linspace(0,pi/2,Naz/4 + 1);

% fault vectors in NW+Z:
% fault plane
fnorm = [-sind(dip)*sind(str);
         -sind(dip)*cosd(str);
          cosd(dip)];
% auxiliary plane normal: slip vector
slipv = [ cosd(rak)*cosd(str) + sind(rak)*cosd(dip)*sind(str);
         -cosd(rak)*sind(str) + sind(rak)*cosd(dip)*cosd(str);
          sind(rak)*sind(dip)];
  
str1 = d2r(str);
dip1 = d2r(dip);

dip2 = acos(sind(rak)*sind(dip));
% str2 = d2r(str - acosd(-1/(tand(dip)*tan(dip2))));
str2 = mod(atan2(slipv(1),slipv(2)),2*pi);
if dip2 > pi/2
    dip2 = pi - dip2;
    str2 = pi + str2;
    str2 = mod(str2,2*pi);
end
if ang_diff(str1,str2) < pi/2
    str2 = pi + str2;
    str2 = mod(str2,2*pi);
end

% points for fault plane
D1   = asin(sin(faz).*sin(dip1));
phi1 = atan(tan(faz).*cos(dip1));
x1_ = sin(phi1).*(cos(D1/2) - sin(D1/2))*rmax; x1_ = [x1_  fliplr(x1_(1:end-1))];
y1_ = cos(phi1).*(cos(D1/2) - sin(D1/2))*rmax; y1_ = [y1_ -fliplr(y1_(1:end-1))];
x1 = cos(str1)*x1_ + sin(str1)*y1_;
y1 = -sin(str1)*x1_ + cos(str1)*y1_;
% points for aux. plane
D2   = asin(sin(faz).*sin(dip2));
phi2 = atan(tan(faz).*cos(dip2));
x2_ = sin(phi2).*(cos(D2/2) - sin(D2/2))*rmax; x2_ = [x2_  fliplr(x2_(1:end-1))];
y2_ = cos(phi2).*(cos(D2/2) - sin(D2/2))*rmax; y2_ = [y2_ -fliplr(y2_(1:end-1))];
x2 = cos(str2)*x2_ + sin(str2)*y2_;
y2 = -sin(str2)*x2_ + cos(str2)*y2_;      

% edge
x3 = rmax.*cos(az);
y3 = rmax.*sin(az);

figure(fign);
subplot(subp(1),subp(2),subp(3));
pcolor(X,Y,-Z(:,:,psvsh))
caxis([-0.5 0.5])
colormap(redbluecmap)
shading interp
axis equal
box off; axis off
text(0,4.5,sprintf('%s radiation pattern',upper(char(phases(psvsh)))),'Fontsize',14,'HorizontalAlignment','center')
hold on
plot(x1,y1,'k','LineWidth',2)
plot(x2,y2,'k','LineWidth',2)
plot(x3,y3,'k','LineWidth',3)

end

