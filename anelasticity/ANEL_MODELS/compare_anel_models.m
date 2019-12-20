%% Script to compare a suite of models
clear
addpath('~/Documents/MATLAB/geoff/VpVsQ/')
addpath('JF10Fit')
%% Figure 1 
% replicate figure 1 from Abers et al., 2014

d = [5e-6, 5e-3];
f = [0.01,1];
P = [0.1,2.5];
T_C = [1000:20:1500]'; T = T_C+273;

figure(1);clf
for iff = 1:2
for ip = 1:2

    for it = 1:length(T)
    [ ~,~,~,q_takei(it)] = takei2017 ( T(it),d(ip),P(ip),2*pi*f(iff),1 );
    [ ~,~,~,q_mcc(it)] = mccarthy2011( T(it),d(ip),P(ip),2*pi*f(iff),1 );
    [ ~,~,~,q_jf10(it)] = JF2010     ( T(it),d(ip),P(ip),2*pi*f(iff),1 );
    [ ~,~,~,q_fj05(it)] = FJ2005     ( T(it),d(ip),P(ip),2*pi*f(iff),1 );
    [ ~,~,~,q_pm13(it)] = PM2013     ( T(it),d(ip),P(ip),2*pi*f(iff),1 );
    end

    subplot(2,2,iff+2*(ip-1))
    hold on
    
    plot(T_C,1./q_takei)
    plot(T_C,1./q_mcc,'r')
    plot(T_C,1./q_jf10,'b')
    plot(T_C,1./q_fj05,'m')
    plot(T_C,1./q_pm13,'g')
end
end
% 
subplot(2,2,1), ylim([0 30])
subplot(2,2,2), ylim([0 100])
subplot(2,2,3), ylim([0 200])
subplot(2,2,4), ylim([0 200])

%% Figure 2
%  replicate Figure 10 in Takei AREPS
T_C = [500:50:1350]'; T = T_C+273;
Z = [50, 75];
for iz = 1:length(Z)
for it = 1:length(T)

[ ~,~,~,q(it,iz),gg(it,iz)] = takei2017 ( T(it),1e-3,Z(iz)/29.94,2*pi/100,1 );
end
end
Vs = sqrt(gg/3350)/1000;

figure(2); clf
subplot(311),semilogy(T_C,q(:,1),'-k',T_C,q(:,2),'--k'), set(gca,'xlim',[1000,1400],'ylim',[1e-4 0.1])
subplot(3,1,2:3),semilogy(T_C,Vs(:,1),'-b',T_C,Vs(:,2),'-m')
