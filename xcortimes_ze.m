function [dcor,dcstd,dvcstd,acor,Gmat,avec]=xcortimes_ze(data,dt,varargin)
% [dcor,dcstd,dvcstd,acor,Gmat,avec]=xcortimes_ze(data,dt,...)
%         something sort of like the vandecar and crosson cross correlation
%
%  REQUIRED INPUTS
%   data = (npt x nsta) array of data
%   dt = 1/samplerate
% 
%  OPTIONAL INPUTS: (as 'name',val pairs)
%   'lagmax' = max allowed lag (in seconds, default=1) 
%   'ifplot' = option for plots (T/F, default=false)
%   't0'     = time of first sample relative to 0 lag (i.e. -ive if before)
%   'ifwt'   = option to weight least squares by xcor of each pair (T/F,
%              default - false)
%   'gclim'  = optional maximum distance between stations to use data; must
%              go along with lat,lon details for all stations. 
%   'scaleopt'= 'unbiased' (by default) or 'biased', or 'none'
%
% dcstd = Geoff's version of the stds
% dvcstd = stds from Vandecar & Crosson eqn 8
% acor = mean correlation coefficient of each trace (with stack)
% 
% NEGATIVE DCOR MEANS EARLIER ARRIVAL, POSITIVE MEANS LATER 
% 
% Z. Eilon, 2018 (adapted from script by G. Abers)

% defaults:
lagmax = 1; % in seconds
ifplot = false; % don't plot
scaleopt = 'unbiased';
t0 = 0; % in seconds
ifwt = false; % unweighted inversion
latlon = [];
gclim = [];


% parse varargin
narginchk(2,inf);
iv = 1;
while iv < nargin-2
    switch varargin{iv}
        case 'lagmax' 
            lagmax = varargin{iv+1};
        case 'ifplot'
            ifplot = varargin{iv+1};
        case 't0'
            t0 = varargin{iv+1};
        case 'ifwt'
            ifwt = varargin{iv+1};
        case 'gclim'
            gclim = varargin{iv+1};           
        case 'latlon'
            latlon = varargin{iv+1};           
        case 'scaleopt'
            scaleopt = varargin{iv+1};   
    end
    iv = iv+2;
end
    

    

nsta=size(data,2);
npt=size(data,1);
nmat=nsta*(nsta-1)/2; % total number of combinations of statinos

% Gmat=sparse(nmat,nsta);
% sparse G matrix making
si = zeros(nmat,1);
sj = zeros(nmat,1);
ss = zeros(nmat,1);

avec=zeros(nmat,1);
wvec=ones(nmat,1);
lmax=round(lagmax./dt);   % Max lag allowed, should set in fcn

% groups of stations
Ngrp = 1;
grps = cell({});

kk=0;
%if (iplot==1) disp('building lag matrix...'); end
for ii=1:(nsta-1) % for each station...
    for jj=(ii+1):nsta % loop over combos with all stations after it
        %  check if stations too far apart?
        if ~isempty(gclim)
            if distance(latlon(ii,1),latlon(ii,2),latlon(jj,1),latlon(jj,2)) > gclim
                continue 
            end
        end
        % do xcorr
        kk=kk+1;
        [c,lags]=xcorr(data(:,ii),data(:,jj),lmax,scaleopt);
        pk=find(c==max(c));
        avec(kk) = lags(pk(1));
        wvec(kk) = max(c);
        si(2*kk+[-1,0]) = kk*[1 1]; % i indices of G
        sj(2*kk+[-1,0]) = [ii jj];  % j indices of G
        ss(2*kk+[-1,0]) = [1 -1];  % j indices of G
        
        % assign to groups
        if kk == 1, grps{1}(1)=ii; end
        for ig = 1:Ngrp
            if ismember(ii,grps{ig})
                % ii is in group ig; assign jj to that group and get out
                grps{ig} = unique([grps{ig};jj]);
                break 
            end
            if ismember(jj,grps{ig})
                % jj is in group ig; assign ii to that group and get out
                grps{ig} = unique([grps{ig};ii]);
                break 
            end
            % if at this point, ii and jj are in no group, make a new group
            Ngrp = Ngrp+1;
            grps{Ngrp} = [ii;jj];
        end
        
%         Gmat(kk,ii) = 1;
%         Gmat(kk,jj) = -1;
    end
end
% adjust vectors for non-pairs
avec = avec(1:kk);
wvec = wvec(1:kk);
si = si(1:2*kk); sj = sj(1:2*kk); ss = ss(1:2*kk);
Gmat = sparse(si,sj,ss,kk,nsta);

% make sure groups have no overlap
grped = 0;
while grped == 0
    regrp = 0;
    for ig1 = 1:Ngrp
        if regrp; break; end
        for ig2 = ig1+1:Ngrp
            if regrp; break; end
            if any(intersect(grps{ig1},grps{ig2})) % if groups have overlapping members
                grps{ig1} = unique([grps{ig1}(:);grps{ig2}(:)]); % put all members in first group
                grps(ig2) = []; % kill second group
                Ngrp = Ngrp-1; % 1 fewer groups
                regrp = 1; % break cycle to restart with new groupings
            end
        end
    end
    if regrp; continue; end
    grped = 1; % can only get here if no regrouping because all groups have unique stations
end
% use only largest group
Ngrpmax = 0;
for ig = 1:Ngrp
    if length(grps{ig})>Ngrpmax
        igrpmax = ig; 
        Ngrpmax = length(grps{ig});
    end
end  

% get out if only two stations
if Ngrpmax==2
    dcor=nan(nsta,1);
    dcstd=nan(nsta,1);
    dvcstd=nan(nsta,1);
    acor=nan(nsta,1);
    return
end

% find stations with no info
noinfstas = find(sum(abs(Gmat)) == 0);
% include stations from all but largest group
for ig = 1:Ngrp
    if ig~=igrpmax
        noinfstas = union(noinfstas(:),grps{ig}(:));
    end
end   

% now for any stas with no info, weak norm damp to zero
Nnostas = length(noinfstas);
H = sparse(1:Nnostas,noinfstas,ones(Nnostas,1),Nnostas,nsta);
h = zeros(Nnostas,1);

% fix mean to zero - per group
D = zeros(Ngrp,nsta);
d = zeros(Ngrp,1);
for ig = 1:Ngrp
    D(ig,grps{ig}) = 1; 
end
D = sparse(D);

% assemble
Gmat = [Gmat;D;H];
avec   = [avec;d;h];
wvec = [wvec;ones(Ngrp+Nnostas,1)];


% 
% 
% % rule out stations too far apart?
% if ~isempty(gclim)
%     kk = 0;
%     dist = zeros(nmat+1,1);
%     for ii=1:(nsta-1) % for each station...
%     for jj=(ii+1):nsta % loop over combos with all stations after it
%         kk=kk+1;
%         dist(kk) = distance(latlon(ii,1),latlon(ii,2),latlon(jj,1),latlon(jj,2));
%     end
%     end
%     % hugely down-weight distant stations
%     wvec(dist>gclim) =  wvec(dist>gclim)*1e-9;
%     % now for any stas with no info, norm damp to zero
%     noinfstas = find(abs(sum(diag(wvec)*Gmat)-1) < 1e-8);
%     Nnostas = length(noinfstas);
%     H = sparse(1:Nnostas,noinfstas,ones(Nnostas,1),Nnostas,nsta);
%     h = zeros(Nnostas,1);
%     Gmat = [Gmat;H];
%     avec   = [avec;h];
%     wvec = [wvec;ones(Nnostas,1)];
% end

% weight by max xcor?
if ifwt
    wmat = diag(wvec); 
else
    wmat = eye(length(wvec));
end
    
    

%if (iplot==1) disp('Inverting lags...'); end
dcor = (Gmat'*wmat*Gmat)\Gmat'*wmat*avec;

resid=avec(1:kk) - Gmat(1:kk,:)*dcor;
dcstd=sqrt(diag(inv(Gmat'*Gmat)).*var(resid));

% Find standard dev of residuals using Vandecar & Crosson equation 8
% R = zeros(nsta) ;
% R(tril(true(nsta),-1)) = resid(1:(end-1-Nnostas));
% R = R - R';
% V = R.^2; 
% dvcstd = sqrt(sum(V,2)./(nsta-2));
for is = 1:nsta
    dvcstd(is,1) = sqrt(sum(resid(Gmat(1:kk,is)~=0).^2)./(sum(Gmat(1:kk,is)~=0)-2));
end

% Find acor
stak = zeros(npt,1);
tstak = [0:npt-1]';
for is=1:nsta
    stak=stak+interp1(tstak-dcor(is),data(:,is),tstak);
end
stak=stak./nsta;
jmx=npt-ceil(max(dcor));
jmn=ceil(abs(min(dcor)))+1;
ida=find(isfinite(stak) & (1:npt)'<jmx & (1:npt)'>jmn );
acor = zeros(nsta,1);
for is=1:nsta
    acor(is)=sum(stak(ida).*data(ida+round(dcor(is)),is))./(std(stak(ida),1).*std(data(ida+round(dcor(is)),is),1))./length(ida);
end

dcor = dcor .* dt;
dcstd=full(dcstd).* dt;

%% Nan stations with no info
dcor(noinfstas) = nan;
dcstd(noinfstas) = nan;
acor(noinfstas) = nan;
dvcstd(noinfstas) = nan;

%% plotting
if (ifplot==1) 
    figure(46)
    clf
    %disp('...finishing up');
tt=dt.*(0:(npt-1))+t0; %edited zje: was tt0=dt.*(0:(npt-1))-pretime;
for ii=1:nsta
    subplot(211);
    plot(tt,data(:,ii),'Linewidth',1.5);
    title('Pre-alignment')
    hold on;
    subplot(212);
    plot(tt-dcor(ii),data(:,ii),'Linewidth',1.5);
    title('Post-alignment')
    hold on;
end
end
return
