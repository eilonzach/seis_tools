% Manipulating moment tensors
% plots fault plane solutions

str = 303.45; % strike of FAULT plane in degrees
dip = 86.17; % dip of fault in degrees
rak = -177.89; % rake of slip (+90 is pure thrust, 0 is L-lateral)
Mo = 1.0530e+25; % seismic moment

% Coordinate system:
coords = 'ZSE'; % common options: NWZ, NE-Z, ZSE, ENZ

fprintf('Strike: %.1f%s  Dip: %.1f%s  Rake: %.1f%s\n\n',...
    str,char(176),dip,char(176),rak,char(176));

%% Part 1: Calculate moment tensor elements in mode coordinate system
% (r, theta, phi) = (up, south, east).
% use M_ij = s_i*n_j + s_j*n_i
% find slip vector in mode coord system 
s = [1; 0; 0]; % slip in 1 direction
n = [0; 1; 0]; % normal in 2 direction
M = s*n' + n*s';

% consider as a set of rotation matrices transforming from the fault coord
% system which has slip vector (1,0,0) and fault plane normal (0,1,0)
% thus transforming to the strike-dip-rake system is a series of negative
% rotations: 
% 1) +rak about [0,1,0] % R-lateral conventionally +ive - changes sign
% 2) -(dip about [1,0,0]
% 3) -str about [0,0,1] 

R1 = [cosd(rak) 0 sind(rak); 0 1 0; -sind(rak) 0 cosd(rak)];
R2 = [1 0 0; 0 -sind(dip) -cosd(dip); 0 cosd(dip) -sind(dip)];
R3 = [cosd(str) sind(str) 0; -sind(str) cosd(str) 0; 0 0 1];
R = R3*R2*R1; % puts it into N W Z

if strcmp(coords,'NE-Z') % then to go from  [N W Z] to [N E -Z]
R4 = [1 0 0; 0 -1 0; 0 0 -1];
elseif strcmp(coords,'ZSE') % or to go from [N W Z] to [Z S E]
R4 = [0 0 1; -1 0 0; 0 -1 0];
elseif strcmp(coords,'NWZ') % or to keep as [N W Z]
R4 = eye(3);
elseif strcmp(coords,'ENZ') % or to go from [N W Z] to [E N Z]
R4 = [0 -1 0; 1 0 0; 0 0 1];
elseif strcmp(coords,'ENZ') % or to go from [N W Z] to [N -Z W]
R4 = [1 0 0; 0 0 -1; 0 -1 0];
end

R_ = R4*R;
sr = R_*s;
nr = R_*n;
fprintf('In %s coordinate system, ',coords)
Mr = Mo*R_*M*R_'

%% Part 2: Calculate moment tensor elements in ENZ coordinate system
R4 = [0 -1 0; 1 0 0; 0 0 1];
R_ = R4*R;
sr = R_*s;
nr = R_*n;
fprintf('In ENZ coordinate system, ',coords)
Mr = Mo*R_*M*R_'

%% Part 3: Calculate and plot radiation pattern
% now do everything in NE-Z coord system
R_ = [1 0 0; 0 -1 0; 0 0 -1]*R;
Mr = Mo*R_*M*R_';
[U,D] = eig(Mr);
D=diag(D);
is1 = find(D==max(D));
is3 = find(D==min(D));
is2 = setdiff([1:3],[is1,is3]);

tr1 = r2d(atan2(U(2,is1),U(1,is1))); pl1 = 90 - acosd(U(3,is1));
if pl1<0, tr1 = mod(tr1+180,360); pl1 = -pl1; end
tr2 = r2d(atan2(U(2,is2),U(1,is2))); pl2 = 90 - acosd(U(3,is2));
if pl2<0, tr2 = mod(tr2+180,360); pl2 = -pl2; end
tr3 = r2d(atan2(U(2,is3),U(1,is3))); pl3 = 90 - acosd(U(3,is3));
if pl3<0, tr3 = mod(tr3+180,360); pl3 = -pl3; end

fprintf('eig1: %.f trend: %.1f%s plunge: %.1f%s\n',D(is1),tr1,char(176),pl1,char(176));
fprintf('eig2: %.f trend: %.1f%s plunge: %.1f%s\n',D(is2),tr2,char(176),pl2,char(176));
fprintf('eig3: %.f trend: %.1f%s plunge: %.1f%s\n',D(is3),tr3,char(176),pl3,char(176));


% figure(10); clf; hold on
% Nr = 8;
% r = linspace(0,4,Nr+1)';
% for ir = 1:length(r)
%     inc = 2*asin(sind(45)*r(ir)/max(r)); %incidence angle
%     az = linspace(0,2*pi,6*(ir-1) + 1)';
%     if ir==1, az = [0;0]; end
%     for iz = 1:length(az)-1
%         p = [sin(inc)*cos(az(iz));...
%              sin(inc)*sin(az(iz));...
%              cos(inc)];
%         u = p'*Mr*p; % instantaneous motion (scalar)
%         h = polar(az(iz),r(ir),'o');
%         set(h,'MarkerSize',10*abs(u))
%         if u>0
%             set(h,'MarkerFaceColor','b')
%         end
%     end    
% end
% polar(d2r(tr1),4*cosd(pl1),'ro')
% polar(d2r(tr2),4*cosd(pl2),'ro')
% polar(d2r(tr3),4*cosd(pl3),'ro')
% polar(0,2,'ko')
% hold off

Nr = 50;
Naz = 360;
r = linspace(0,4,Nr)'; % each column a different r
inc = 2*asin(sind(45)*r/max(r));
% for polar plot, r = r_max * (sin(inc/2)/sin(45))
%                   = r_max * (cos(dip/2) - sin(dip/2))
az = linspace(0,2*pi,Naz+1)'; % each row a different az
X = sin(az)*r';
Y = cos(az)*r';
Z = zeros(size(X));
for ir = 1:length(r)
for iz = 1:length(az)
    p = [sin(inc(ir))*cos(az(iz));...
         sin(inc(ir))*sin(az(iz));...
         cos(inc(ir))];
    Z(iz,ir) = p'*Mr*p; % instantaneous motion (scalar)
end    
end
figure(11); clf
pcolor(X,Y,-Z)
caxis(Mo*[-0.01 0.01])
colormap('gray')
shading interp
axis equal
title('Focal Mechanism in lower hemisphere projection')

hold on
plot(4*cos([0:0.01:2*pi]),4*sin([0:0.01:2*pi]),'-k') %border
text(4*(cosd(pl1/2)-sind(pl1/2))*sind(tr1),4*(cosd(pl1/2)-sind(pl1/2))*cosd(tr1),...
    ['\fontsize{16}\bf{\color{blue}T}'],'BackgroundColor','w','Margin',7,'HorizontalAlignment','center')
text(4*(cosd(pl2/2)-sind(pl2/2))*sind(tr2),4*(cosd(pl2/2)-sind(pl2/2))*cosd(tr2),...
    ['\fontsize{16}\bf{\color{green}N}'],'BackgroundColor',[0.5 0.5 0.5],'Margin',7,'HorizontalAlignment','center')
text(4*(cosd(pl3/2)-sind(pl3/2))*sind(tr3),4*(cosd(pl3/2)-sind(pl3/2))*cosd(tr3),...
    ['\fontsize{16}\bf{\color{red}P}'],'BackgroundColor','k','Margin',7,'HorizontalAlignment','center')
hold off

%% Part 4: Moment tensor looking along strike
% Starting off in [N W Z], looking vertically
%  rotate back by strike around the vertical axis
R_ = [cosd(str) -sind(str) 0; sind(str) cosd(str) 0; 0 0 1]*R;
% and switch to N-ZW
R_ = [1 0 0; 0 0 -1; 0 1 0]*R_;
Mr = Mo*R_*M*R_';

sr = [cosd(rak); -sind(rak)*sind(dip); sind(rak)*cosd(dip)];
nr = [0; -cosd(dip); -sind(dip)];
fprintf('In strike, down, left coordinate system, ')
Mr = sr*nr' + nr*sr'

%% Part 5: Calculate and plot back-hemisphere projection in this orientation
Nr = 50;
Naz = 360;
r = linspace(0,4,Nr)'; % each column a different r
inc = 2*asin(sind(45)*r/max(r));
% for polar plot, r = r_max * (sin(inc/2)/sin(45))
%                   = r_max * (cos(dip/2) - sin(dip/2))
az = linspace(0,2*pi,Naz+1)'; % each row a different az
X = sin(az)*r';
Y = cos(az)*r';
Z = zeros(size(X));
for ir = 1:length(r)
for iz = 1:length(az)
    p = [cos(inc(ir));...
         -sin(inc(ir))*cos(az(iz));...
         -sin(inc(ir))*sin(az(iz))];
    Z(iz,ir) = p'*Mr*p; % instantaneous motion (scalar)
end    
end
figure(12); clf
pcolor(X,Y,-Z)
caxis([-0.01 0.01])
colormap('gray')
shading interp
hold on
plot(4*cos([0:0.01:2*pi]),4*sin([0:0.01:2*pi]),'-k') %border
hold off
axis equal
title('FPS looking along strike')
