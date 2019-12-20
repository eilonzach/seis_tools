clear; close all
addpath ANEL_MODELS/

figure(1),clf,set(gcf,'pos',[440 87 668 711])

T_K = linspace(600,1326,1000) + 273; 
depth_km = 50;
P_GPa = depth_km/29.94;
d_mm = 1;
phi = 0;
period_sec = 100;
for j=1:length(T_K)
    [j1,j2] = takei2017(T_K(j),d_mm/1000,P_GPa,2*pi/period_sec,1,[],phi);
    J50(j,2) = j1 + 1i*j2;
end

rho  = 3300;
V    = 1./sqrt(rho*real(J50));
Qinv = imag(J50)./real(J50);
    
subplot(3,1,2:3); plot(T_K-273,V/1000,'b'); ylabel('V_S, km/sec'); xlabel('T ^\circ C');  grid on; hold on;
plot((1326)*[1 1],[4 5],'--b');
subplot(3,1,1); semilogy(T_K-273,Qinv,'-k'); ylabel('Q^{-1}'); xlabel('T ^\circ C');  grid on; hold on;


T_K = linspace(850,1351,1000) + 273; 
depth_km = 75;
P_GPa = depth_km/29.9401;
d_mm = 1;
phi = 0;
period_sec = 100;
for j=1:length(T_K)
    [j1,j2] = takei2017(T_K(j),d_mm/1000,P_GPa,2*pi/period_sec,1,[],phi);
    J75(j,2) = j1 + 1i*j2;
end

rho  = 3300;
V    = 1./sqrt(rho*real(J75));
Qinv = imag(J75)./real(J75);
    
subplot(3,1,2:3); plot(T_K-273,V/1000,'m'); axis square; axis([400 1500 4.15 4.65]);
plot((1351)*[1 1],[4 5],'--m');
subplot(3,1,1); semilogy(T_K-273,Qinv,'--k'); axis([1000 1400 10^-4 10^-1]);
