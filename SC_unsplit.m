function [phiSC,dTSC, Ru, Tu, Emap, inipolest,phiEV,dTEV,Lmap] = SC_unsplit( datN, datE, samprate, plotopt )
% [phiSC,dTSC, Ru, Tu, Emap, inipolest,phiEV,dTEV,Lmap] = SC_unsplit( datN, datE, samprate, plotopt )
% do Silver and Chan unsplitting
%
% INPUTS:
%     datN      = time series of data on N component
%     datE      = time series of data on E component
%     samprate  = sample rate of data (per second)
%     plotopt   = option to plot (1) or not (0)
% 
% OUTPUTS:
%     phiSC     = estimated fast direction of anisotropy (Min Energy method)
%     dTSC      = estimated splitting time (Min Energy method)
%     Ru        = time series of unsplit data on R component
%     Tu        = time series of unsplit data on T component
%     Emap      = energy map from grid search (Min Energy method)
%     inipolest = estimated initial polarisation of wave
%  - if nargout < 7 ignore this
%     phiEV     = estimated fast direction of anisotropy (Eigenvalue method)
%     dTEV      = estimated splitting time (Eigenvalue method)     
%     Lmap      = eigenvalue map from grid search (Eigenvalue method)

if nargin < 4
    plotopt = 0;
end

ddT  = 0.1;
dphi = 2;

%% Measure initial polarisation following Vidale 1986
[ inipolest ] = pol_vidale_simple2D( datN,datE );

%% rotate data to R and T
[ datRT ] = zne2zrt( [datN,datE], inipolest );
%% Grid search to unsplit
phis = -90:dphi:90;
dTs  =   0:ddT :4;
Emap = zeros(length(phis),length(dTs)); % Minimum energy method
Lmap = zeros(length(phis),length(dTs)); % Eigenvalue method

for ip = 1:length(phis) 
for it = 1:length(dTs)
    [ Ru,Tu ] = unsplit( datRT(:,1), datRT(:,2), phis(ip), dTs(it), samprate, inipolest );
    Cu = [Ru'*Ru, Ru'*Tu; Tu'*Ru, Tu'*Tu];
    Emap(ip,it) = Cu(2,2);
    if nargout>6 % if also doing EV method
        [~,Du] = eigs(Cu);
        Lmap(ip,it) = Du(2,2);
    end
end
end

%Minimum energy method
[~,it,ip] = mingrid(Emap);
phiSC = phis(ip);
dTSC = dTs(it);

%Eigenvalue method
if nargout>6 
    [~,it2,ip2] = mingrid(Lmap);
    phiEV = phis(ip2);
    dTEV = dTs(it2);
end

[ Ru,Tu ] = unsplit( datRT(:,1), datRT(:,2), phiSC, dTSC, samprate, inipolest );

if plotopt
    figure(1); clf
    subplot(3,1,1)
    plot([datN,datE]), title('datN (blue) datE (green)')
    subplot(3,1,2)
    plot(datRT), title('datR (blue) datT (green)')
    subplot(3,1,3)
    plot([Ru,Tu]), title('datR-unsplit (blue)  datT-unsplit (green)')
    
    figure(2)
    pcolor(Emap)
    set(gca,'Xtick',[1:0.5/ddT:length(dTs)],'XtickLabel',dTs(1:0.5/ddT:end))
    set(gca,'Ytick',[1:30/dphi:length(phis)],'YtickLabel',phis(1:30/dphi:end))
    set(gca,'FontSize',12)
    hold on
    plot(it,ip,'wo','Linewidth',3)
%     plot(it2,ip2,'ko','Linewidth',3)
end

end

