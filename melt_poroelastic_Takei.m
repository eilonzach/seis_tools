function [ Vpu_m, Vsu_m ] = melt_poroelastic_Takei( Vpu,Vsu,phi,A,Km )
% [ Vpu_m, Vsu_m ] = melt_poroelastic_Takei( Vpu,Vsu,phi,[A=1.6], [Km=30e9] )
% 
%  THIS FUNCTION IS ENTIRELY COPIED FROM BEN + CHRIS' VBR FUNCTION:
%       el_ModUnrlx_MELT_f.m
%  ALL CREDIT AND THANKS (+BLAME) TO BEN HOLTZMAN AND CHRIS HAVLIN
%       
%  ZE 2016
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poro-elastic effect of melt. 
% The effect of melt fraction is calculated using Takei's 2002 Appendix A 
% fit to the isotropic solutions of her generalized solution in her 1995 
% paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUTS:
%  Vpu = unrelaxed Vp (m/s)
%  Vsu = unrelaxed Vs (m/s)
%  phi = melt fraction
%  A   = wetting angle factor (1:2.3, Yoshino)
%  Km  = bulk modulus of the melt (Pa)

%% prelims
if nargin < 4 || isempty(A)
    A = 1.6 ; % 1:2.3 depending upon the wetting angle (see Yoshino). 
end
if nargin < 5 || isempty(Km)
    Km = 30e9; % melt bulk modulus [Pa], Takei 2002, Table 2 
end

    

rho = 3300; % this value does not matter for this - it cancels!
Gu = rho*Vsu.^2;
Mu = rho*Vpu.^2;
Ku = Mu - (4/3)*Gu;
nu = (Mu - 2*Gu)./(2*Mu - 2*Gu); %Poisson's ratio

%% calculate effective moduli and standard deviations
for ii = 1:length(Gu)
  [Gu_eff,Gamma_G(ii,1)]=melt_shear_moduli(Gu(ii),phi(ii),A,nu(ii)) ;
  [Ku_eff,Gamma_K(ii,1)]=melt_bulk_moduli(Ku(ii),phi(ii),A,Km,nu(ii));    
end

%% calculate effective velocities and standard deviations
   [Vpu_m,Vsu_m] = Vp_Vs_calc(phi,Gu,nu,Gamma_G,Gamma_K,rho,Km); 
  

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates Vp and Vs, accouting for poro-elastic effects. Reduces to 
% pure phase calculation when phi = 0. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vp,Vs] = Vp_Vs_calc(phi,Gu,nu,Gamma_G,Gamma_K,rho,K_m)

  delta = 1e-20; % small number to avoid division by 0.
  
% unrelaxed bulk mod  
  nu_fac = 2/3*(1+nu)./(1-2*nu); 
  Ku = Gu.*nu_fac; 
  
% the poroelastic factors  
  K1 = (1-Gamma_K).^2; 
  K2 = 1 - phi - Gamma_K + phi.*Ku./K_m; % equals 0 if phi = 0; 
  
% effective bulk and shear modulus  
  bulk_mod = Ku .* (Gamma_K + K1./(K2 + delta)); 
  shear_mod = Gu .* Gamma_G; 

% calculate Vp, Vs  
  [Vp,Vs] = VpVs_from_moduli(bulk_mod,shear_mod,rho);
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% melt_bulk_moduli
% calculates the bulk moduli, accounting for poro-elastic effect of melt. 
% following Takei [citation]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kb_eff,Gamma_k]=melt_bulk_moduli(k,phi,A,Km,nu) 

  
% calculate the contiguity as a function of the melt fraction
  Psi =  1-A.*sqrt(phi) ;    

% ==========================================================
% exponents n and m in Takei, 2002 (functions of contiguity)
% coefficients:

  a = zeros(3,4); % this will look like the transpose of what Takei has in table A1: 
  a(1,:) = [1.8625 0.52594 -4.8397 0 ];
  a(2,:) = [4.5001 -6.1551 -4.3634 0 ];
  a(3,:) = [-5.6512 6.9159 29.595 -58.96];
  
  sz_a = size(a);
  a_vec=zeros(1,sz_a(1));   
  for i = 1:sz_a(1)
      a_vec(i) = sum(a(i,:).*nu.^(0:1:3)); 
  end  
  
% eqns A5, A6, Takei 2002
  n_k = a_vec(1).*Psi + a_vec(2).*(1-Psi) + a_vec(3).*(1-Psi).^1.5 ;

% ==================================
% normalized skeleton properties as functions of Psi

  Gamma_k = (1-phi).*(1-(1-Psi).^n_k).*ones(size(k));   
  k_sk_prime = 1-(1-Psi).^n_k ; % eqn A3
  K_sk = k_sk_prime.*k ; % eqn A3
  
% =====================================================
% the effective laws with melt:
  delta =1e-20 ; % small number to avoid dividing by zero at phi=0; 
  top = (1-K_sk./k).^2 ;
  bot = (1-phi-K_sk./k + phi.*k./Km) + delta.*(phi==0) ; 
  Kb_eff_prime = K_sk./k + top./bot ;
  Kb_eff = Kb_eff_prime.*k;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% melt_shear_moduli
% calculates the shear moduli, accounting for poro-elastic effect of melt. 
% following Takei [citation]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mu_eff,Gamma_Mu]=melt_shear_moduli(mu,phi,A,nu) 

  
% calculate the contiguity as a function of the melt fraction
  Psi =  1-A.*sqrt(phi) ;  
  
% FIX THIS !!:w check this ! relation between mu and K !!
% calculate Vp/Vs as the outcome of this too !!! 

% ==========================================================
% exponents n and m in Takei, 2002 (functions of contiguity)
% coefficients:

  b = zeros(3);
  b(1,:) = [1.6122 0.13572 0] ;
  b(2,:) = [4.5869 3.6086 0] ;
  b(3,:) = [-7.5395 -4.8676 -4.3182] ;  
  
  sz_b = size(b);
  b_vec=zeros(1,sz_b(1)); 
  for i = 1:sz_b(1)      
      b_vec(i) = sum(b(i,:).*nu.^(0:1:2)); 
  end  

% eqns A6, Takei 2002
  n_mu = b_vec(1).*Psi + b_vec(2).*(1-Psi) + b_vec(3).*(1-Psi).^2 ;

% ==================================
% normalized skeleton properties as functions of Psi
  Gamma_Mu=(1-phi).*(1-(1-Psi).^n_mu).*ones(size(mu)); 
  mu_sk_prime = 1-(1-Psi).^n_mu  ; % eqn A4
  Mu_sk = (1-phi).*mu_sk_prime.*mu ; % eqn A4

% =====================================================
% the effective laws with melt:

  Mu_eff = Mu_sk ; % effective shear modulus
  
  
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VpVs_from_moduli
% calculates the Vp and Vs from the bulk and shear moduli and rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vp,Vs] = VpVs_from_moduli(bulk_mod,shear_mod,rho)

Vp = sqrt((bulk_mod + (4/3)*shear_mod)./rho);
Vs = sqrt(shear_mod./rho);

end

