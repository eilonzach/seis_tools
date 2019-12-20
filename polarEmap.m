function [ phiest,dTest ] = polarEmap( Emap,R,th,a,dof,figconf,plotconf )
% [ phiest,dTest ] = polarEmap( Emap,R,th,a,dof,figconf,plotconf )
% 
%   This function creates a polar map of error showing errors over a space
%   of points defined by a polar coordinate system. The error within this
%   space is contoured.
% 
% INPUTS: need first 2
%   Emap - error map, XxY with rows for each azimuth and columns for each magnitude.
%   R    - vector of magnitudes corresponding to each column. If this
%           vector is just one element long, assume this is the maximum
%           magnitude (Rmax) and mags are spaced equally from 0 to Rmax
%   th   - vector of azimuths corresponding to each row. If this argument
%           is ommitted then assume they are equally spaced from 90 to
%           -89.999 (i.e. just more than -90 ; e.g. -88:2:90 or -89:1:90)
%   a    - confidence level to plot on
%   dof  - degrees of freedom for this measurement
%   figconf - 3-element vector as follows:
%                (1) - figure number
%                (2) - vertical no. of subplots
%                (3) - horizontal no. of subplots
%   plotconf - n-element vector describing all the subplot spaces occupied
%               by the plot
% 
% OUTPUTS
%  phiest = [minphi avphi maxphi] using confidence interval a
%  dTest  = [mindT  avdT  maxdT ] using confidence interval a

% ie = 8;
% ir = 1;
% R = 4;
% nargin=2;
% Emap = eq(ie).results(ir).EEVmap;
% eq(ie).results(ir).dtEV
% eq(ie).results(ir).phiEV


x = size(Emap,1);
y = size(Emap,2);
if nargin < 6
    figconf=[1,1,2];
    plotconf = 2;
end
if nargin < 5 || any(isnan(dof))==1
    dof=5;
end
if nargin < 4 || any(isnan(a))==1
    a=0.05;
end
if nargin < 3 || any(isnan(th))==1
    dth = 180/x;
%     th = flipud([-90+dth:dth:90]');
    th = ([-90:dth:90]');
else
    dth = unique(diff(th));
end

if length(R)==1
    Rmax = R;
    dR = Rmax/(y-1);
    R = [0:dR:Rmax]';
else
    Rmax = max(R);
    dR = mean(unique(diff(R)));
end

% extend to wrap around to -90
% Emap(end+1,:) = Emap(1,:);
% th(end+1) = -90;

% find minimum
[Emin,ix,iy] = mingrid(Emap);

% mapping -90 to 90 range over whole circle: stretch by doubling azimuths
th = 2*th;

[TH,RR] = meshgrid(d2r(th),R);
[X,Y] = pol2cart(pi/2 - TH,RR); % change so 0 at north, clockwise...
% various contours
clev = mingrid(Emap):(maxgrid(Emap)-mingrid(Emap))/50:maxgrid(Emap); % evenly spaced contours
alev = mingrid(Emap)*(1 + (2/(dof-2))*finv(1-a,dof-2,2)); % 95% contour
olev = alev*(1:0.5:5); % multiples of 95% contour

figure(figconf(1));

if nargin<6 % only plot standard error map if this is a one-off plot
    clf
subplot(1,2,1), cla, hold on
% surf(X,Y,Emap')
%% Standard error map - cartesian
pcolor(R,th./2,Emap); shading interp
plot(R(ix),th(iy)/2,'ow','MarkerSize',8,'MarkerFaceColor','r')
contour(R,th./2,Emap,alev,'w','Linewidth',1);
xlabel('dt (sec)');
ylabel('phi (\circ)');
set(gca,'YTick',-90:30:+90); 
set(gca,'XTick',0:1:4);
end

%% Polar error map
subplot(figconf(2),figconf(3),plotconf)
% polarPcolor(R,th,Emap)
p = polar([0 2*pi], [0 Rmax]);  
delete(p);
gca;

hold on
% hpa = polaraxes; hold on
[C,h] = contourf(X,Y,Emap',clev,'linestyle','none'); 
% [C,h] = contour(X,Y,Emap',olev,'k');
[C,h] = contour(X,Y,Emap',[alev alev],'k','Linewidth',1.5);
h = polar(d2r(90 - th(iy)),R(ix),'ok');
set(h,'MarkerSize',10,'MarkerFaceColor','r')
hold off

% view(90,-90)

%% relabel angles
oldlab = [0:30:330]';
newlab = mod(((90 - oldlab)./2) + 90,180)-90;
for ii = 1:length(oldlab)
    if oldlab(ii) == 270
        newstr = '-90+';
    else
        newstr = num2str(newlab(ii));
        if newlab(ii)>0
            newstr = strcat('+',newstr);
        end
    end
    set(findall(gcf,'String',num2str(oldlab(ii))),'String',newstr,'fontsize',15)
end

%% relabel times
dtlab = [1:4]';
for ii = 1:length(dtlab)
	hold on
    set(findall(gcf,'String',['  ',num2str(dtlab(ii))]),'String','')
    if ii<4
        plot(dtlab(ii)*sin([0:0.01:2*pi]),dtlab(ii)*cos([0:0.01:2*pi]),':k','linewidth',1.5)
        text((dtlab(ii)+0.2)*sin(-pi/6),(dtlab(ii)+0.2)*cos(-pi/6),[num2str(dtlab(ii)),' s'],...
            'fontsize',15,'horizontalalignment','right');
    else 
        plot(Rmax*sin([0:0.01:2*pi]),Rmax*cos([0:0.01:2*pi]),'k','linewidth',2)
    end
    hold off
end


%% Pretty plot
hold on
% plot([0 0 nan -Rmax Rmax],[-Rmax Rmax nan 0 0],'--','linewidth',1.2,'color',[0.5 0.5 0.5])
%% Get minima
[~, minx, miny] = mingridc(Emap,(1 + (2/(dof-2))*finv(1-a,dof-2,2)));
   dTest = (minx-1)*dR;
   phiest = miny*dth - 90;

end


