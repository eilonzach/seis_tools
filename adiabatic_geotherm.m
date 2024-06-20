function [T_C] = adiabatic_geotherm(z_km,T_C_pot,k_opt)
% [T_C] = adiabatic_geotherm(z_km,T_pot[1350],k_opt['vary'])
% 
% function to produce temperatures on an adiabatic geotherm, T_C (defalt
% Celsius) as a function of depth, z_km. Can alternatively input a single
% depth value, in which case output will be for all depths from zero
% (surface) to this point. Can also select k_opt, which is a flag to use a
% constant ('const') value for the dT/dz_ad or a varying value ('vary') to
% account for temperature and depth-dependent thermal conductivity. 'vary '
% is the default.
arguments
    z_km (:,1)
    T_C_pot (1,1) = 1350;
    k_opt = 'vary'
end

zz = 0:2800;

% make vector of T gradients
dT_dz_ad = nan(size(zz));
switch k_opt
    case 'vary'
        % use approx values from Katsura et al. 2010, https://doi.org/10.1016/j.pepi.2010.07.001
        % these values are loosely grabbed from their Fig3
        um = zz <= 410;
        dT_dz_ad(um) = interp1([0 410]',[0.64 0.36]',zz(um));
        utz = (zz > 410) & (zz <= 520);
        dT_dz_ad(utz) = interp1([410 520]',[0.37 0.35]',zz(utz));
        ltz = (zz > 520) & (zz <= 670);
        dT_dz_ad(ltz) = interp1([520 670]',[0.42 0.38]',zz(ltz));
        ulm = (zz > 670) & (zz <= 1500);
        dT_dz_ad(ulm) = interp1([670 1500]',[0.47 0.35]',zz(ulm));        
        llm = (zz > 1500) & (zz <= 2800);
        dT_dz_ad(llm) = interp1([1500 2800]',[0.35 0.24]',zz(llm));  
    case 'const'
        dT_dz_ad(:) = 0.4;
end

% calculate T increments down adiabat
dT = zeros(size(zz));
dT(2:end) = midpts(dT_dz_ad).*diff(zz);

TT_C = T_C_pot + cumsum(dT);

T_C = interp1(zz,TT_C,z_km);





end