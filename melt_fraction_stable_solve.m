function [solution] = melt_fraction_stable_solve(T_K,P_GPa,C_h2o_solid,D_h2o,C_co2_solid,D_co2,solidus_str)
% Using equations in appendix of Blatter et al. (2022), solve for
% thermodynamically stable melt fraction given temperature
%
% In this version, assume melt free and dry solidus below 
% ~330 km depth (P > 11.25 GPa)
%
% INPUTS
% T_K : Temperature in Kelvins [N x 1]
% P_GPa : Pressure in GPa [N x 1]
% C_h2o : ppm H2O contained in solid matrix (i.e., bulk content) [N x 1]
% D_h2o : H2O partition coefficient [1 x 1]
% C_co2 : ppm CO2 contained in solid matrix (i.e., bulk content) [N x 1]
% D_co2 : CO2 partition coefficient [1 x 1]
% solidus_str : Type of solidus solution to use | 'katz' or 'hirschmann'
%
% J.B. Russell 2023

% Cf_CO2_saturation = 38; % (wt %), saturation value used in Blatter et al. (2022)
Cf_CO2_saturation = 37; % (wt %), saturation value to be compatible with Dasgupta et al. (2007)

% Range over which to evaluate F
F_test = [0:0.0001:1]';

solution.phi = zeros(length(P_GPa),1);
solution.Cf_H2O = zeros(length(P_GPa),1);
solution.Cf_CO2 = zeros(length(P_GPa),1);
solution.Cs_H2O = zeros(length(P_GPa),1);
solution.Cs_CO2 = zeros(length(P_GPa),1);
solution.Cs_H2O_0 = zeros(length(P_GPa),1);
solution.Cs_CO2_0 = zeros(length(P_GPa),1);
solution.Cs_H2O_ppm = zeros(length(P_GPa),1);
solution.Cs_CO2_ppm = zeros(length(P_GPa),1);
solution.Tsolidus_K = zeros(length(P_GPa),1);

% Initilize with dry solidus
[dry] = SoLiquidus(P_GPa*1e9,0,0,solidus_str);
solution.Tsolidus_K = dry.Tsol(:) + 273;

for iz = 1:length(P_GPa)
    
    if mod(iz,10000)==0
        disp([num2str(iz),'/',num2str(length(P_GPa))])
    end
    
    if P_GPa(iz) > 11.25 % outside range of Hirschmann (2010) fit in Fig 1
        continue
    end

    Cs_H2O=C_h2o_solid(iz)*1e-4; % ppm to wt%
    Cs_CO2=C_co2_solid(iz)*1e-4; % ppm to wt%

    % calculate volatile fractions in melt
    Cf_H2O = Cs_H2O ./ (D_h2o + F_test * (1-D_h2o));
    Cf_CO2 = Cs_CO2 ./ (D_co2+ F_test * (1-D_co2));
    Cf_CO2(Cf_CO2>Cf_CO2_saturation) = Cf_CO2_saturation;
    Cf_CO2(isnan(Cf_CO2)) = 0;

    P_Pa = P_GPa(iz)*1e9; % [Pa]
    [Solidus] = SoLiquidus(P_Pa,Cf_H2O,Cf_CO2,solidus_str);
    Tsolidus_K = Solidus.Tsol(:) + 273;

    % Hirschmann et al. (2009, 2010) PEPI
    dT_dF = -40*P_GPa(iz) + 450;
%     F = (T_K(iz)-Tsolidus_K) ./ dT_dF; % thermodynamic melt fraction (weight %)
    h = F_test(:) - (T_K(iz)-Tsolidus_K(:)) ./ dT_dF;
    
    % Check whether pressure is higher than Hirschmann's (2010) fit in Fig 1, resulting in negative dT_dF
    Hirschmann10_undefined = 0;
    if dT_dF<0 % corresponds to P > 11.25 GPa
        Hirschmann10_undefined = 1;
    end 
    
%     % Katz et al (2003) Eq. 19
%     h = F_test(:) - ((T_K(iz)-Tsolidus_K(:))./(Solidus.Tlherz_dry(:)-Solidus.Tsol_dry(:))).^1.5;
    
    % Find value of F where h is zero
    F_solution = interp1(h,F_test,0);
    if isnan(F_solution) || Hirschmann10_undefined
        % No melt predicted
        F_solution = 0;
    end
    
    % Save out solution values
    solution.phi(iz) = F_solution;
    solution.Cf_H2O(iz) = Cs_H2O ./ (D_h2o + F_solution * (1-D_h2o));
    solution.Cf_CO2(iz) = Cs_CO2 ./ (D_co2 + F_solution * (1-D_co2));
    if solution.Cf_CO2(iz) > Cf_CO2_saturation
        solution.Cf_CO2(iz) = Cf_CO2_saturation;
    elseif isnan(solution.Cf_CO2(iz))
        solution.Cf_CO2(iz) = 0;
    end
    solution.Cs_H2O(iz) = D_h2o*solution.Cf_H2O(iz);
    solution.Cs_CO2(iz) = D_co2*solution.Cf_CO2(iz);
    solution.Cs_H2O_0(iz) = Cs_H2O;
    solution.Cs_CO2_0(iz) = Cs_CO2;
    solution.Cs_H2O_ppm(iz) = solution.Cs_H2O(iz)*1e4;
    solution.Cs_CO2_ppm(iz) = solution.Cs_CO2(iz)*1e4;
    [Solidus] = SoLiquidus(P_Pa,solution.Cf_H2O(iz),solution.Cf_CO2(iz),solidus_str);
    solution.Tsolidus_K(iz) = Solidus.Tsol(:) + 273;
    
    
    if 0
        figure(999);
        if iz == 1
            clf; 
        end
        hold on;
        plot(F_test,h);
        plot([1 1]*F_solution,[min(h) max(h)],'--k');
        pause;
    end
end

end

