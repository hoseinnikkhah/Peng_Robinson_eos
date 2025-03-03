%% PR-KM Model for Ceftriaxone Sodium Solubility in SC-CO2
% This script implements the Peng-Robinson equation of state with 
% Kwak-Mansoori mixing rules to model the solubility of Ceftriaxone sodium
% in supercritical CO2 according to the paper

clear;
clc;
close all;

%% Constants
R = 8.314; % Gas constant [J/(mol·K)]

%% Input Parameters
% Temperature range
T = [308 318 328 338]; % [K]

% Pressure range
P = 1:1:300; % [bar] from 1 to 300 with step of 1

% Binary interaction parameters from the paper (Figure 4)
k_ij = [0.103 0.117 0.126 0.134]; % Binary interaction parameters for each temperature

% Critical properties of Ceftriaxone sodium
T_c_drug = 1149.5; % Critical temperature [K]
P_c_drug = 10.57;  % Critical pressure [bar]
omega_drug = 2.0964; % Acentric factor

% Critical properties of CO2
T_c_CO2 = 304.18; % Critical temperature [K]
P_c_CO2 = 73.8;   % Critical pressure [bar]
omega_CO2 = 0.225; % Acentric factor

% Parameters for calculating sublimation pressure
T_ref = 304.18;     % Reference temperature in K
P_ref = 1e-6;       % Reference pressure in bar
delta_H_vap = 56910;% Enthalpy of vaporization in J/mol (56.91 kJ/mol)
R_gas = 8.314;      % Gas constant in J/(mol·K)

% Calculate solid molar volume (MW/density)
MW_drug = 598.5;    % Molecular weight of Ceftriaxone sodium in g/mol
density_solid = 1.5;% Assumed density in g/cm³ (typical for pharmaceutical solids)
v_solid = MW_drug/density_solid; % Molar volume in cm³/mol
v_solid_SI = v_solid * 1e-6; % Convert cm³/mol to m³/mol

%% Initialize arrays to store results
num_pressures = length(P);
solubility_values = zeros(4, num_pressures); % Store solubilities for 4 temperatures
density_CO2 = zeros(4, num_pressures);      % Store CO2 densities

%% Calculate PR-KM parameters for pure components
a_CO2 = zeros(1, 4);
b_CO2 = zeros(1, 4);
d_CO2 = zeros(1, 4);
a_drug = zeros(1, 4);
b_drug = zeros(1, 4);
d_drug = zeros(1, 4);

for i = 1:4
    % Calculate EOS parameters for CO2
    [a_CO2(i), b_CO2(i), d_CO2(i)] = calculate_PR_params(T_c_CO2, P_c_CO2, T(i), omega_CO2);
    
    % Calculate EOS parameters for Ceftriaxone sodium
    [a_drug(i), b_drug(i), d_drug(i)] = calculate_PR_params(T_c_drug, P_c_drug, T(i), omega_drug);
end

%% Main calculation loop for each temperature and pressure
for temp_idx = 1:4 % Loop through temperatures
    % Get NIST density data for CO2 or calculate it (you can implement this part)
    % For now, I'll use a simple approximation based on ideal gas law corrected by Z
    
    for p_idx = 1:num_pressures % Loop through pressures
        current_P = P(p_idx);
        current_T = T(temp_idx);
        
        % Calculate sublimation pressure using Clausius-Clapeyron equation
        P_sublimation = P_ref * exp((delta_H_vap/R_gas)*(1/T_ref - 1/current_T));
        
        % Convert pressure from bar to Pa for calculation
        P_current_Pa = current_P * 1e5;
        P_sublimation_Pa = P_sublimation * 1e5;
        
        % Calculate mixture parameters for PR-KM EOS
        % Define cross-interaction parameters
        a_ij = sqrt(a_CO2(temp_idx)*a_drug(temp_idx))*(1-k_ij(temp_idx));
        b_ij = (b_CO2(temp_idx) + b_drug(temp_idx))/2;
        d_ij = (d_CO2(temp_idx) + d_drug(temp_idx))/2;
        
        % Start with initial guess for mole fraction of drug (y_1)
        % Implement iterative method to find the mole fraction that satisfies equation 4
        % For the initial implementation, I'll try a range of possible y_1 values
        
        y_1_values = linspace(1e-10, 1e-5, 100); % Range of potential solubility values to try
        residuals = zeros(size(y_1_values));
        
        for y_idx = 1:length(y_1_values)
            y_1 = y_1_values(y_idx);
            y_2 = 1 - y_1;
            
            % Calculate mixture parameters according to KM mixing rules (equations A.10-A.12)
            b_m = (y_1^2 * b_CO2(temp_idx)) + (2*y_1*y_2*b_ij) + (y_2^2 * b_drug(temp_idx));
            d_m = (y_1^2 * d_CO2(temp_idx)) + (2*y_1*y_2*d_ij) + (y_2^2 * d_drug(temp_idx));
            
            c1 = (y_1^2 * a_CO2(temp_idx)^(2/3) * b_CO2(temp_idx)^(1/3));
            c2 = (2*y_1*y_2 * (a_ij^(2/3)) * (b_ij^(1/3)));
            c3 = (y_2^2 * a_drug(temp_idx)^(2/3) * b_drug(temp_idx)^(1/3));
            c4 = (y_1^2 * b_CO2(temp_idx));
            c5 = (2*y_1*y_2*b_ij);
            c6 = (y_2^2 * b_drug(temp_idx));
            
            a_m = (c1 + c2 + c3)^(1.5) / sqrt(c4 + c5 + c6);
            
            % Calculate A, B, D parameters for PR EOS
            [A, B, D] = calculate_ABD(current_T, current_P, a_m, b_m, d_m);
            
            % Solve the PR EOS for compressibility factor Z
            Z = solve_PR_cubic(A, B, D);
            
            % Calculate fugacity coefficient using equation A.14
            sum_y_b_1j = y_1*b_CO2(temp_idx) + y_2*b_ij;
            sum_y_a_1j = y_1*a_CO2(temp_idx) + y_2*a_ij;
            sum_y_d_1j = y_1*d_CO2(temp_idx) + y_2*d_ij;
            
            tau_1 = (a_m + R*current_T*d_m) - (2*sqrt(a_m*d_m*R*current_T)*(1/2 - (sum_y_b_1j/b_m))) + sum_y_a_1j*(1 - sqrt((R*current_T*d_m)/a_m)) + sum_y_d_1j*(R*current_T - sqrt((R*current_T*a_m)/d_m));
            
            % Calculate terms for fugacity coefficient
            term1 = (Z-1)*((2*sum_y_b_1j/b_m) - 1);
            term2 = -log(Z - B);
            term3 = -(tau_1/(sqrt(2)*R*current_T*b_m))*log((Z + (1 + sqrt(2))*B)/(Z + (1 - sqrt(2))*B));
            
            ln_phi = term1 + term2 + term3;
            phi = exp(ln_phi);
            
            % Calculate Poynting factor
            poynting_factor = exp((v_solid_SI * (P_current_Pa - P_sublimation_Pa)) / (R_gas * current_T));
            
            % Calculate right-hand side of equation 4
            rhs = (P_sublimation * poynting_factor) / (current_P * phi);
            
            % Calculate residual (difference between assumed and calculated y_1)
            residuals(y_idx) = abs(y_1 - rhs);
        end
        
        % Find the y_1 value that minimizes the residual
        [~, min_idx] = min(residuals);
        solubility_values(temp_idx, p_idx) = y_1_values(min_idx);
        
        % For plotting purposes, we also need CO2 density at each condition
        % This would typically come from NIST data or another equation of state
        % As a placeholder, I'll estimate it using the ideal gas law with Z correction
        Z_CO2 = 0.8; % Approximate Z for supercritical CO2 - in practice, this would be calculated
        density_CO2(temp_idx, p_idx) = (current_P * 1e5 * 44.01) / (Z_CO2 * R_gas * current_T * 1000); % kg/m³
    end
end

%% Plot solubility vs pressure (like Figure 4)
figure;
set(gcf, 'Position', [100, 100, 800, 600]);

% Plot symbols for experimental data (you would need to add the actual experimental data points)
markers = ['o', 's', '^', 'd']; % Markers for different temperatures
colors = ['k', 'r', 'g', 'b'];  % Colors for different temperatures

% Plot model curves
for i = 1:4
    semilogy(P, solubility_values(i, :)*1e6, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', 2);
    hold on;
end

% Add labels and legend
xlabel('Pressure (bar)', 'FontSize', 12);
ylabel('Drug solubility y \times 10^6', 'FontSize', 12);
legend('T = 308 K, k_{ij} = 0.103, PR EoS/KM', ...
       'T = 318 K, k_{ij} = 0.117, PR EoS/KM', ...
       'T = 328 K, k_{ij} = 0.126, PR EoS/KM', ...
       'T = 338 K, k_{ij} = 0.134, PR EoS/KM', ...
       'Location', 'NorthWest');
title('Solubility of Ceftriaxone sodium in SC-CO2 (PR-KM Model)', 'FontSize', 14);
grid on;
axis([0 300 0.01 10]);

%% Plot solubility vs density (for comparison with Figure 3b)
figure;
set(gcf, 'Position', [100, 100, 800, 600]);

% Plot solubility vs CO2 density
for i = 1:4
    semilogy(density_CO2(i, :), solubility_values(i, :)*1e6, 'LineStyle', '-', 'Color', colors(i), 'LineWidth', 2);
    hold on;
end

% Add labels and legend
xlabel('CO2 Density (kg/m³)', 'FontSize', 12);
ylabel('Drug solubility y \times 10^6', 'FontSize', 12);
legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K', 'Location', 'NorthWest');
title('Solubility of Ceftriaxone sodium vs CO2 Density', 'FontSize', 14);
grid on;

%% Helper functions
function [a, b, d] = calculate_PR_params(T_c, P_c, T, omega)
    % Calculate parameters for PR EOS based on critical properties
    R = 8.314; % J/(mol·K)
    
    % Convert P_c from bar to Pa
    P_c_Pa = P_c * 1e5;
    
    % Calculate m parameter
    m = 0.37464 + 1.5422*omega - 0.26992*omega^2;
    
    % Calculate a, b parameters
    a = 0.457235 * (R^2 * T_c^2 / P_c_Pa) * (1 + m * (1 - sqrt(T/T_c)))^2;
    b = 0.077796 * (R * T_c / P_c_Pa);
    
    % Calculate d parameter for Kwak-Mansoori
    d = 0.457235 * (R^2 * T_c^2 / P_c_Pa) * (1 + m * (1 - sqrt(T/T_c)))^2 * (m^2) / (R*T_c);
end

function [A, B, D] = calculate_ABD(T, P, a, b, d)
    % Calculate dimensionless parameters for PR EOS
    R = 8.314; % J/(mol·K)
    
    % Convert P from bar to Pa
    P_Pa = P * 1e5;
    
    A = a * P_Pa / ((R * T)^2);
    B = b * P_Pa / (R * T);
    D = d * P_Pa / (R * T);
end

function Z = solve_PR_cubic(A, B, D)
    % Solve PR cubic equation for Z using Newton-Raphson method
    % Coefficients for cubic equation: Z³ + c₂Z² + c₁Z + c₀ = 0
    
    % Calculate cubic equation coefficients
    c2 = B - 1;
    c1 = D - 3*B^2 - 2*B + A - 2*sqrt(A*D);
    c0 = B^3 + B^2 - A*B - B*D + 2*B*sqrt(A*D);
    
    % Initial guess for Z
    Z = 1.0;
    
    % Newton-Raphson iteration
    max_iter = 100;
    tol = 1e-8;
    
    for i = 1:max_iter
        % Evaluate function and derivative
        f = Z^3 + c2*Z^2 + c1*Z + c0;
        df = 3*Z^2 + 2*c2*Z + c1;
        
        % Update Z
        Z_new = Z - f/df;
        
        % Check convergence
        if abs(Z_new - Z) < tol
            Z = Z_new;
            break;
        end
        
        Z = Z_new;
    end
    
    % Make sure Z > B (physically meaningful)
    if Z <= B
        Z = B + 1e-6;
    end
end
