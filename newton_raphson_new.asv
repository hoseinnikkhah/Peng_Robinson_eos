function newton_raphson_new()
    run('PR_Mansoori.m');
    
    % Constants and properties for Ceftriaxone sodium
    MW_drug = 598.5;                % Molecular weight of Ceftriaxone sodium (g/mol)
    rho_solid = 1.3;                % Approximate density of solid drug (g/cm³) - estimated
    v_solid = (MW_drug/rho_solid)/1000;  % Molar volume of solid (m³/mol), converted to L/mol
    
    % Parameters for sublimation pressure calculation
    Delta_H_vap = 56.91 * 1000;     % Heat of vaporization (J/mol), converted from kJ/mol
    T_ref = 293;                    % Reference temperature (K)
    P_ref = 1e-10;                  % Reference pressure at T_ref (bar) - very low for pharmaceutical compounds
    
    % Settings
    T_fixed = 308;                  % Fixed temperature for the test
    P_range = 1:1:300;              % Pressure range (bar)
    
    % Arrays to store results
    calculated_y = zeros(size(P_range));
    Z_values = zeros(size(P_range));
    phi_values = zeros(size(P_range));
    
    for idx = 1:length(P_range)
        P_fixed = P_range(idx);
        
        % Calculate sublimation pressure at T_fixed using Clausius-Clapeyron equation
        P_sub = P_ref * exp((Delta_H_vap/R) * ((1/T_ref) - (1/T_fixed)));
        
        % Initial guess for y_1
        y_1 = 1e-6;  % Start with a very small value typical for drug solubility
        
        % Calculate mixture properties based on the guess
        y_2 = 1 - y_1;

        b_m = (y_1^2 * b_ij_3D(1,1,1)) + (2*y_1*y_2*b_ij_3D(1,2,1)) + (y_2^2 * b_ij_3D(2,2,1));
        d_m = (y_1^2 * d_ij_3D(1,1,1)) + (2*y_1*y_2*d_ij_3D(1,2,1)) + (y_2^2 * d_ij_3D(2,2,1));

        c1 = (y_1^2 * a_ij_3D(1,1,1)^(2/3)*b_ij_3D(1,1,1)^(1/3));
        c2 = (2*y_1*y_2 * (a_ij_3D(1,2,1)^(2/3)) * (b_ij_3D(1,2,1)^(1/3)));
        c3 = (y_2^2 * a_ij_3D(2,2,1)^(2/3) * b_ij_3D(2,2,1)^(1/3));
        c4 = (y_1^2 * b_ij_3D(1,1,1));
        c5 = (2*y_1*y_2*b_ij_3D(1,2,1));
        c6 = (y_2^2 * b_ij_3D(2,2,1));

        a_m = (c1 + c2 + c3)^(1.5) / sqrt(c4 + c5 + c6);

        [A, B, D] = correlations_cap(T_fixed, P_fixed, a_m, b_m, d_m);

        sum_y_b_1j = y_1*b_ij_3D(1,1,1) + y_2*b_ij_3D(1,2,1);
        sum_y_a_1j = y_1*a_ij_3D(1,1,1) + y_2*a_ij_3D(1,2,1);
        sum_y_d_1j = y_1*d_ij_3D(1,1,1) + y_2*d_ij_3D(1,2,1);

        tau_1 = (a_m + R*T_fixed*d_m) - (2*sqrt(a_m*d_m*R*T_fixed)*(1/2 - (sum_y_b_1j/b_m))) + sum_y_a_1j*(1 - sqrt((R*T_fixed*d_m)/a_m)) + sum_y_d_1j*(R*T_fixed - sqrt((R*T_fixed*a_m)/d_m));
      
        % Try multiple initial guesses to find a physically meaningful root
        initial_guesses = [0.2, 0.5, 0.8]; % Try liquid-like, mid-range, and gas-like
        best_root = NaN;
        best_phi = NaN;

        for guess_idx = 1:length(initial_guesses)
            initial_guess = initial_guesses(guess_idx);
            max_iterations = 100;
            tolerance = 1e-10;
            
            [root, iterations, converged] = newton_raphson_cubic(initial_guess, max_iterations, tolerance, A, B, D);
            
            % Check if this is a valid root
            if converged && isfinite(root) && root > B
                best_root = root;
                
                % Calculate phi using the calculated root Z
                coeff1 = 2*(y_1*b_ij_3D(1,1,1) + y_2*b_ij_3D(1,2,1));
                
                % Protect against denominator in coeff2 being too small
                denom = root + (1 - sqrt(2))*B;
                if abs(denom) < 1e-10
                    denom = sign(denom) * 1e-10;
                    if denom == 0
                        denom = 1e-10;
                    end
                end
                
                coeff2 = (root + (1 + sqrt(2))*B)/denom;
                
                % Calculate ln(phi)
                term1 = (root-1)*((coeff1/b_m) - 1);
                term2 = -log(root - B);
                term3 = -(tau_1/(sqrt(2)*R*T_fixed*b_m))*log(coeff2);
                
                ln_phi = term1 + term2 + term3;
                
                if ln_phi > -700 && ln_phi < 700  % Check for reasonable ln(phi) value
                    best_phi = exp(ln_phi);
                    break;  % We found a good root, no need to try more
                end
            end
        end
        
        % Calculate the solubility using the thermodynamic relationship
        if isfinite(best_root) && isfinite(best_phi) && best_phi > 0
            % Assume phi_saturated = 1 (common assumption for pharmaceutical solids)
            phi_saturated = 1.0;
            
            % Poynting factor
            poynting_factor = exp(v_solid * (P_fixed - P_sub) / (R * T_fixed));
            
            % Calculate solubility
            y_calculated = (phi_saturated * P_sub * poynting_factor) / (P_fixed * best_phi);
            
            % Store results
            calculated_y(idx) = y_calculated;
            Z_values(idx) = best_root;
            phi_values(idx) = best_phi;
            
            fprintf('P = %d bar: Z = %.6f, phi = %.6e, y_calculated = %.6e\n', ...
                    P_fixed, best_root, best_phi, y_calculated);
        else
            fprintf('Warning: Could not find valid solution at P = %d bar\n', P_fixed);
            calculated_y(idx) = NaN;
            Z_values(idx) = NaN;
            phi_values(idx) = NaN;
        end
    end
    
    % Filter out any NaN values for plotting
    valid_indices = ~isnan(calculated_y);
    
    % Create a figure to match Figure 4 from the paper
    figure;
    plot(P_range(valid_indices), calculated_y(valid_indices) * 1e6, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Pressure (bar)');
    ylabel('Drug Solubility (y × 10^6)');
    title('Solubility of Ceftriaxone Sodium in SC-CO2 at T = 308 K');
    
    % Save results to CSV
    results_table = table(P_range', Z_values', phi_values', calculated_y', ...
                          'VariableNames', {'Pressure', 'Z', 'Phi', 'Solubility'});
    writetable(results_table, 'solubility_results.csv');
    fprintf('Results saved to solubility_results.csv\n');
end

function [root, iterations, convergence] = newton_raphson_cubic(initial_guess, max_iter, tolerance, A, B, D)
    z = initial_guess;
    iterations = 0;
    convergence = false;
    
    sqrt_AD = sqrt(A*D);
    
    coef_z2 = B - 1;
    
    coef_z1 = D - 3*B^2 - 2*B + A - 2*sqrt_AD;
    
    coef_z0 = B^3 + B^2 - A*B - B*D + 2*B*sqrt_AD;
    
    for i = 1:max_iter
        f_z = z^3 + coef_z2*z^2 + coef_z1*z + coef_z0;
        
        f_prime_z = 3*z^2 + 2*coef_z2*z + coef_z1;
        
        if abs(f_prime_z) < 1e-10
            f_prime_z = sign(f_prime_z) * 1e-10;
            if f_prime_z == 0
                f_prime_z = 1e-10;
            end
        end
        
        z_next = z - f_z / f_prime_z;
        
        error = abs(z_next - z);
        
        z = z_next;
        iterations = i;
        
        if error < tolerance
            convergence = true;
            break;
        end
    end
    
    root = z;
end