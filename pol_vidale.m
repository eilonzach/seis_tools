function [ P ] = pol_vidale( dat,tt0, avwind,bazi )
%[ P ] = pol_vidale( dat,tt0, avwind,bazi )
% function to estimate the polarisation of a data series on 2/3 component
% seismogram
%
% INPUTS
% dat:      time series of windowed, filtered (decimated) data, 3/2
%           components, in columns [RR,TT(,ZZ)]
% tt0:      epochal times of each of the datapoints (REMEMBER DECIMATION)
% avwind:   number of points in covariance averaging window
% OUTPUTS
% P.time:   start time of the window
% P.strike: direction of polarisation

%% setup
dt=mean(diff(tt0));
polD=size(dat,2);
if polD==2 
	dat=[dat,zeros(size(dat,1),1)]; % pad out to a Nx3 with zeros ==> ZZ=0
elseif polD < 2 || polD > 3
    error('Need data for 2 horizontal or all 3 components')
end

%% Measure polarisation angle
% make structure for polarisation    
n=size(dat,1)-avwind+1;
P=struct('time',{tt0(1)+ dt*([0:1:n-1]'+(avwind-1)/2)},...
        'dip',zeros(n,1),'phiRT',zeros(n,1),'phiNE',zeros(n,1),'L0',zeros(n,1),...
        'Pe',zeros(n,1),'Ps',zeros(n,1),'Pp',zeros(n,1),'Dphi',zeros(n,1));
wind=hanning(avwind);
for j=1:n  % loop over averaging windows

RR=dat(j:j+avwind-1,1).*wind;
TT=dat(j:j+avwind-1,2).*wind;
ZZ=dat(j:j+avwind-1,3).*wind;
% Measure incoming polarisation angle - VIDALE 1986
uu = hilbert(RR);
vv = hilbert(TT);
ww = hilbert(ZZ);
%N.B. uu,vv,ww, are column vectors

C=[uu'*conj(uu) uu'*conj(vv) uu'*conj(ww);...
   vv'*conj(uu) vv'*conj(vv) vv'*conj(ww);...
   ww'*conj(uu) ww'*conj(vv) ww'*conj(ww)];
C=C(1:polD,1:polD); % only use 2x2 top left part for 2D

[V,D] = eigs(C); %eigenstructure in order of abs(eigenvalue)
D=diag(real(D)); %C is Hermitian and therefore has real, positive eigenvalues
ei1=1; %find(D==max(D));
ei3=polD; %find(D==min(D));
if polD==3; 
ei2=2; %find(D~=max(D) && D~=min(D));
end;

P.L0(j)=ei1;
x0=V(1,ei1);
y0=V(2,ei1);
if polD==3, z0=V(3,ei1); else z0=0; end
%search over angles to maximise X
X=zeros(180,1);
for alpha=1:180
    cisa = cosd(alpha) + 1i*sind(alpha);
    X(alpha) = sqrt(real(x0*cisa)^2 + real(y0*cisa)^2 + real(z0*cisa)^2);
end
[X,alpha]=max(X); %N.B. reassignment of X and alpha to their max-X values
% rotate the vector (x0,y0,z0) by alpha in the x-y plane
rot =  [ cosd(alpha) sind(alpha); -sind(alpha) cosd(alpha)];
x0=rot(1,:)*V(1:2,ei1);
y0=rot(2,:)*V(1:2,ei1);
% z0=rot(3,:)*V(:,ei1);
%% Estimate elliptical component of polarisation
% Pe is 1 for circularly polarized motion and 0 for linearly polarised
% motion
Pe = sqrt(1-X^2)/X;
P.Pe(j)=Pe;
% can find the sense of rotation by comparing the signs of the real and
% imaginary parts of the rotated eigenvector
%% Calculate strike of direction of max. polarisation
phi=r2d(-atan2(real(y0),real(x0))); % NOT SURE ABOUT THE MINUS SIGN
P.phiRT(j)=phi;
P.phiNE(j)=mod(bazi+phi, 360);

if imag(y0)/real(y0) > 1 && imag(x0)/real(x0) > 1
    disp('Note, Im(y0) and Im(x0) > Re(y0) and Re(x0)');
end
%% Calculate dip of direction of max. polarisation
if polD==3
del=r2d(atan2(real(z0),sqrt(real(x0)^2 + real(y0)^2)));
P.dip(j)=del;
end
%% Estimate strength of polarisation
%Ps is near 1 if the signal is completely polarised, i.e. only one major
%component of polarisation - 0 otherwise
if polD==3
Ps = 1 - (D(ei2)+D(ei3))/D(ei1);
else
Ps = 1 - D(ei3)/D(ei1);
end
P.Ps(j)=Ps;
%% Estimate degree of planar polarisation
% Pp is 1 if intermediate component of polarisation is much larger than the
% smallest component, 0 if D2 and D3 are comparable
if polD==3
Pp = 1 - D(ei3)/D(ei2);
P.Pp(j)=Pp;
end

end % loop on averaging windows

end % function
