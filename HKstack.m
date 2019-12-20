function [Amap,H_vec,K_vec] = HKstack(RF,tt,rayp,phase_wts,Vs_av,H_vec,K_vec)
% Amap = HKstack(RF,tt,Vs_av,H_grid,K_grid)

%% grids of H and K
if nargin < 6 || isempty(H_vec)
    H_vec = [20:1:65]';
end
H_vec = H_vec(:); % make column

if nargin < 7 || isempty(K_vec)
    K_vec = [1.6:0.01:1.9];
end
K_vec = K_vec(:)'; % make row

Vs_av = Vs_av.*ones(size(K_vec));


%% compute predicted arr
Vp_av = K_vec.*Vs_av;
kks = sqrt(Vs_av.^-2 - rayp.^2);
kkp = sqrt(Vp_av.^-2 - rayp.^2);
% using equations 2-4 from Zhu and Kanamori 2000
t_ps = H_vec*(kks - kkp); % using outer product
t_ppps = H_vec*(kks + kkp); % using outer product
t_ppss = H_vec*(2*kks); % using outer product

% sum weighted contributions from each phase type
phase_wts = phase_wts/sum(phase_wts);
Amap =  phase_wts(1).*interp1(tt,RF,t_ps) ...
      + phase_wts(2).*interp1(tt,RF,t_ppps) ...
      - phase_wts(3).*interp1(tt,RF,t_ppss,[],0); % negative phase!


end

