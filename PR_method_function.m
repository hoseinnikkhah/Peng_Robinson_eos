function PR_method_function(T_input, P_input)
    run('PR_Mansoori.m');
    
    initial_z_values = 0:10^-5:1;
    num_initial_values = length(initial_z_values);

    roots = zeros(size(initial_z_values));
    phi_values = zeros(size(initial_z_values));
    tau_values = zeros(size(initial_z_values));
    convergence_flags = zeros(size(initial_z_values));
    iteration_counts = zeros(size(initial_z_values));
    calculated_y1_values = zeros(size(initial_z_values));
    relative_differences = zeros(size(initial_z_values));

    T_fixed = T_input;
    P_fixed = P_input;
    
    % Parameters for calculating sublimation pressure
    T_ref = 304.18;        % Reference temperature in K
    P_ref = 1e-6;       % Reference pressure in bar
    delta_H_vap = 56910;% Enthalpy of vaporization in J/mol (56.91 kJ/mol)
    R_gas = 8.314;      % Gas constant in J/(mol·K)
    
    % Calculate sublimation pressure using Clausius-Clapeyron equation
    P_sublimation = P_ref * exp((delta_H_vap/R_gas)*(1/T_ref - 1/T_fixed));
    
    % Calculate solid molar volume (MW/density)
    MW_drug = 598.5;    % Molecular weight of Ceftriaxone sodium in g/mol
    density_solid = 1.5;% Assumed density in g/cm³ (typical for pharmaceutical solids)
    v_solid = MW_drug/density_solid; % Molar volume in cm³/mol
    
    % Convert to m³/mol for calculation
    v_solid_SI = v_solid * 1e-6; % Convert cm³/mol to m³/mol
    
    % Convert bar to Pa for calculation
    P_fixed_Pa = P_fixed * 1e5;
    P_sublimation_Pa = P_sublimation * 1e5;

    for i = 1:num_initial_values
        y_1 = initial_z_values(i);
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
      
        tau_values(i) = tau_1;

        % Try multiple initial guesses to ensure we find a physically meaningful root
        initial_guesses = [0.2, 0.5, 0.8]; % Try liquid-like, mid-range, and gas-like
        best_root = NaN;
        best_iterations = 0;
        best_converged = false;

        for guess_idx = 1:length(initial_guesses)
            initial_guess = initial_guesses(guess_idx);
            max_iterations = 100;
            tolerance = 1e-10;
            
            [root, iterations, converged] = newton_raphson_cubic(initial_guess, max_iterations, tolerance, A, B, D);
            
            % Check if this is a valid root and store the best one
            if converged && isfinite(root)
                if root > B  % Physically meaningful root
                    best_root = root;
                    best_iterations = iterations;
                    best_converged = converged;
                    break;  % We found a good root, no need to try more
                elseif isnan(best_root)  % If we don't have any root yet
                    best_root = root;
                    best_iterations = iterations;
                    best_converged = converged;
                end
            end
        end
        
        % Always use the best root found, even if it's not ideal
        if isnan(best_root)
            roots(i) = NaN;
            phi_values(i) = NaN;
            calculated_y1_values(i) = NaN;
            relative_differences(i) = NaN;
            convergence_flags(i) = 0;
            iteration_counts(i) = max_iterations;
        else
            root = best_root;
            iterations = best_iterations;
            converged = best_converged;
            
            % If root is less than or equal to B, adjust it slightly
            if root <= B
                root = B + 1e-6;  % Ensure we are slightly above B for log(Z-B)
            end
            
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
            
            % Calculate ln(phi) first to check for overflow/underflow
            term1 = (root-1)*((coeff1/b_m) - 1);
            term2 = -log(root - B);
            term3 = -(tau_1/(sqrt(2)*R*T_fixed*b_m))*log(coeff2);
            
            ln_phi = term1 + term2 + term3;
            
            % Handle extreme ln(phi) values to prevent overflow/underflow
            if ln_phi < -700
                phi = 0;
                calculated_y1_values(i) = NaN;
                relative_differences(i) = NaN;
            elseif ln_phi > 700
                phi = Inf;
                calculated_y1_values(i) = NaN;
                relative_differences(i) = NaN;
            else
                phi = exp(ln_phi);
                
                % Calculate Poynting factor
                poynting_factor = exp((v_solid_SI * (P_fixed_Pa - P_sublimation_Pa)) / (R_gas * T_fixed));
                
                % Calculate new y_1
                new_y_1 = (1 * P_sublimation * poynting_factor) / (P_fixed * phi);
                calculated_y1_values(i) = new_y_1;
                
                % Calculate relative difference
                relative_differences(i) = abs(initial_z_values(i) - new_y_1) / new_y_1;
            end
            
            roots(i) = root;
            phi_values(i) = phi;
            convergence_flags(i) = converged;
            iteration_counts(i) = iterations;
        end
    end
    
    % Filter out any NaN values for comparison
    valid_comparison_indices = ~isnan(calculated_y1_values) & ~isnan(roots);
    
    % Save all results to CSV
    results_table = table(initial_z_values', roots', phi_values', tau_values', calculated_y1_values', relative_differences', convergence_flags', iteration_counts', ...
                          'VariableNames', {'InitialMoleFraction', 'Z', 'Phi', 'Tau', 'CalculatedMoleFraction', 'RelativeDifference', 'Converged', 'Iterations'});
    writetable(results_table, 'fugacity_results.csv');
    
    % Create and save the new comparison table with only valid comparisons
    pressure_column = P_fixed * ones(sum(valid_comparison_indices), 1);
    temperature_column = T_fixed * ones(sum(valid_comparison_indices), 1);
    
    comparison_table = table(...
        calculated_y1_values(valid_comparison_indices)', ...
        initial_z_values(valid_comparison_indices)', ...
        relative_differences(valid_comparison_indices)' * 100, ... % Convert to percentage
        pressure_column, ...
        temperature_column, ...
        'VariableNames', {'CalculatedMoleFraction', 'InitialMoleFraction', 'RelativeDifferencePercent', 'Pressure', 'Temperature'});
    
    writetable(comparison_table, 'mole_fraction_comparison.csv');
end
%'
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

function result = conditional(condition, true_value, false_value)
    if condition
        result = true_value;
    else
        result = false_value;
    end
end