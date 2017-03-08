function [dcor,dcstd,dvcstd,acor]=xcortimes(dtrn,dt,pretime,lagmax,ifplot,scaleopt)
% [dcor, dcstd, dvcstd, acor]=xcortimes(dtrn, dt, pretime, lagmax,ifplot,scaleopt)
%         something sort of like the vandecar and crosson cross correlation
%
%  dtrn = (npt x nsta) array of data
%  dt = samplerate
%  pretime = time of first sample relative to 0 lag
%  lagmax = max allowed lag in seconds to test
%  ifplot=1 to do graphics, 0 to just return lags
%  scaleopt = 'ubiased' (by default) or 'biased', or 'none'
%
% dcstd = Geoff's version of the stds
% dvcstd = stds from Vandecar & Crosson eqn 8
% acor = mean correlation coefficient of each trace (with stack)
% 
% NEGATIVE DCOR MEANS EARLIER ARRIVAL, POSITIVE MEANS LATER 
% 
% edits ZE 08/2013
% edit ZE 01/2017 - added [scaleopt='unbiased'] to xcorr routine - CAN
% CAUSE ERRORS IF NOT!!
if nargin < 3 || isempty(pretime)
    pretime=0;
end
if nargin < 4 
    lagmax=[];
end
if nargin < 5 || isempty(ifplot)
    ifplot=0;
end
if nargin < 6 || isempty(scaleopt)
    scaleopt = 'unbiased';
end
    

nsta=size(dtrn,2);
npt=size(dtrn,1);
nmat=nsta*(nsta-1)/2; % total number of combinations of statinos
amat=zeros(nmat+1,nsta);
avec=zeros(nmat+1,1);
lmax=round(lagmax./dt);   % Max lag allowed, should set in fcn
kk=0;
%if (iplot==1) disp('building lag matrix...'); end
for ii=1:(nsta-1) % for each station...
    for jj=(ii+1):nsta % loop over combos with all stations after it
        kk=kk+1;
        [c,lags]=xcorr(dtrn(:,ii),dtrn(:,jj),lmax,scaleopt);
        pk=find(c==max(c));
        avec(kk) = lags(pk(1));
        amat(kk,ii) = 1;
        amat(kk,jj) = -1;
    end
end
amat(kk+1,:)=ones(1,nsta); % damping
%if (iplot==1) disp('Inverting lags...'); end
dcor = amat \ avec;
resid=avec - amat*dcor;
dcstd=sqrt(diag(inv(amat'*amat)).*var(resid));

% Find standard dev of residuals using Vandecar & Crosson equation 8
R = zeros(nsta) ;
R(tril(true(nsta),-1)) = resid(1:end-1);
R = R - R';
V = R.^2; 
dvcstd = sqrt(sum(V,2)./(nsta-2));

% Find acor
stak = zeros(npt,1);
tstak = [0:npt-1]';
for is=1:nsta
    stak=stak+interp1(tstak-dcor(is),dtrn(:,is),tstak);
end
stak=stak./nsta;
jmx=npt-ceil(max(dcor));
jmn=ceil(abs(min(dcor)))+1;
ida=find(isfinite(stak) & (1:npt)'<jmx & (1:npt)'>jmn );
acor = zeros(nsta,1);
for is=1:nsta
    acor(is)=sum(stak(ida).*dtrn(ida+round(dcor(is)),is))./(std(stak(ida),1).*std(dtrn(ida+round(dcor(is)),is),1))./length(ida);
end

dcor = dcor .* dt;
dcstd=dcstd .* dt;

if (ifplot==1) 
    figure(46)
    clf
    %disp('...finishing up');
tt0=dt.*(0:(npt-1))-pretime; %edited zje: was tt0=dt.*(0:(npt-1))-pretime;
for ii=1:nsta
    subplot(211);
    plot(tt0,dtrn(:,ii),'Linewidth',1.5);
    title('Pre-alignment')
    hold on;
    subplot(212);
    plot(tt0-dcor(ii),dtrn(:,ii),'Linewidth',1.5);
    title('Post-alignment')
    hold on;
end
end
return
