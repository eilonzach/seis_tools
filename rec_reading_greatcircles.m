% Quick script to plot great circle paths for Record Readings
% By Stephen A. Veitch LDEO

% It looks best if you save it as a post-script file. 

% You need the mapping toolbox for this to work, but that should be
% included in the full-version of MATLAB (as far as I know).

% The things you should need to modify are in the first section, but, if
% you're really feeling fancy you can play around with the map projection
% (that's on line 58)

%% Modify the variables below to fit desired plot

%Define Generized Map Type:  Matlab needs different inputs for square
   % map projections (ie  Mercator) and for Orhographic Projectiosn
   % the latter require a centre point rather than edge coödinates, and cover
   % ~180? surrounding that point

%define lat/long limits OR center point
latmin=20; % Bottom Latitude
latmax=70;  % Top Latitude

longmin=-120; % Left Longitude
longmax=30; % Right Logitude 

Pline=10; % Grid spacing for latidue lines
Mline=10; % Grid spacing for longitude lines

%define event location
elat=44.799; % Event Latitude
elong=11.192; % Event Longitude

%input station co?rdinates [lat lon] 2-column vector
%current stations available at
%http://www.iris.edu/earthscope/usarray/#transportable_array
TAsta=TAarray; %Name of MATLAB variable containing station locations

% Script will eliminate stations not plotted on the record section
pazmin=60; % Minimum arc dist. included on RS
pazmax=80; % Maximum arc dist. included on RS

% Script, by default, will plot arc distance markers every 10?
% If desired, it will plot finer markers over the area covered 
% in the record section
FINEPLOT=0; % 1 for yes, 0 for no
FINESPACING=1; % Spacing for fine plot

% SHOULDN'T NEED TO MODIFY ANYTHING BELOW HERE
% %
% %
% %
% %

%% 

% Define station locations
slat=TAsta(:,1);
slong=TAsta(:,2);

% plot map
figure('color','w');
hh=axesm('mapproj','eqdconicstd',...
   'maplatlim',[latmin latmax],'maplonlim',[longmin longmax],...
   'MLineLocation',Mline,'PLineLocation',Pline,'Grid','on',...
   'MeridianLabel','on','ParallelLabel','on','MLabelParallel','south');
axis off, gridm on, framem on;
geoshow('landareas.shp')

% Remove stations not included on RS
sdist=acosd(sind(elat).*sind(slat)+cosd(elat)...
   .*cosd(slat).*cosd(slong-elong));
sf1=find((sdist>=pazmin)&(sdist<=pazmax));
slat=slat(sf1);
slong=slong(sf1);

%plot great circle paths
hold on

for i=1:length(slat)
   [la, lo]=gcwaypts(elat,elong,slat(i),slong(i),30);
   geoshow(slat(i),slong(i),'color','k','marker','.');
   hh=geoshow(la,lo,'displaytype','line','color','r');
   set(hh,'LineWidth',1)
end

%Plot event and stations
geoshow(elat,elong,'DisplayType','point',...
   'markeredgecolor','k','markerfacecolor','k','marker','o')
for i=1:length(slat)
   geoshow(slat(i),slong(i),'marker','v',...
       'markerfacecolor','k','markeredgecolor','k');
end

%plot fine distances
if FINEPLOT==1;
   for i=[pazmin:FINESPACING:pazmax];            
   [azlat azlong]=scircle1(elat,elong,i);
   geoshow(azlat,azlong,'color','k','LineStyle','--','LineWidth',.5);
   end
end

%plot coarse distance
azis=[0:1:360];
for i=10:10:170
   [azlat azlong]=scircle1(elat,elong,i);
   geoshow(azlat,azlong,'color','k');
   textm(azlat(25),azlong(25),{i})
end

% Remove variables used
clear latmin latmax longmin longmax Pline Mline elat elong TAsta 
clear pazmin pazmax FINEPLOT FINESPACING slat slong hh i sf1
clear azlat azlong azis la lo sdist