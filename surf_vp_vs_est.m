function [ vp_est,vs_est,E_surf,vp_range,vs_range ] = surf_vp_vs_est( datZRT,tt,rayp,wind,ifplot)
% [ vp_est,vs_est,E_surf,vp_range,vs_range ] = surf_vp_vs_est( datZRT,tt,rayp,windmifplot,ifplot )
% 
%  Function to estimate the Vp and Vs values in the uppermost crust by grid
%  searching through possible values of each parameter. 
% 
% datZRT is in convention: Z+up, R+away, T(not used)
% rayp is in s/km
% wind is a time window [startT endT] that refers to the times in tt

if nargin < 5 || isempty(ifplot)
    ifplot = false;
end

%% prep data
ind = (tt < max(wind)) & (tt >= min(wind)); % find values within window
Zdat = -datZRT(ind,1); % take negative, as Z has to be positive down
Rdat =  datZRT(ind,2);
nsamps = sum(double(ind));

%% rayp
if rayp>1
    rayp = rayp_sdeg2skm(rayp);
end


%% grid search
vp_range = [3.  : 0.025 : 7.5]; Np = length(vp_range);
vs_range = [1.3 : 0.04 : 4.1]; Ns = length(vs_range);

E_surf = nan(Np,Ns);
for is = 1:Ns
for ip = 1:Np

    %condition on Vp/Vs ratio
    vpvs = vp_range(ip)/vs_range(is);
    if vpvs > 2,    continue; end
    if vpvs < 1.40, continue; end
    if isempty(Rdat) || isempty(Zdat), continue; end

    [P,SV] = Rotate_XZ_to_PSV(Rdat,Zdat,vp_range(ip),vs_range(is),rayp);

    % minimise dot product of parent and daughter
    E_surf(ip,is) = abs(P'*SV)/nsamps;
    
%     % minimise energy on daughter
%     E_surf(ip,is) = abs(SV'*SV)/nsamps;

%     % minimise smallest eigenvalue;
%     S = svd([P,SV],0);
%     E_surf(ip,is) = S(2);

end
end
% best guess
[~,x,y] = mingrid(E_surf);
vs_est = vs_range(x);
vp_est = vp_range(y);

if ifplot
    % plot error surface
    figure(4),clf, hold on
    contourf(vs_range,vp_range,E_surf,40,'edgecolor','none'), colorbar
    scatter(vs_est,vp_est,50,'w','filled')

    [Pest,SVest] = Rotate_XZ_to_PSV(Rdat,Zdat,vp_est,vs_est,rayp);

    max1 = max(max(abs([Zdat,Rdat])));
    max2 = max(max(abs([Pest,SVest])));

    figure(2), clf
    subplot(2,2,1), plot(tt(ind),-Zdat,'k'); ylim(max1*[-1 1]); ylabel('Z comp')
    subplot(2,2,3), plot(tt(ind),Rdat,'r');  ylim(max1*[-1 1]); ylabel('R comp')
    subplot(2,2,2), plot(tt(ind),Pest,'k');  ylim(max2*[-1 1]); ylabel('P comp')
    subplot(2,2,4), plot(tt(ind),SVest,'r'); ylim(max2*[-1 1]); ylabel('SV comp')
end

