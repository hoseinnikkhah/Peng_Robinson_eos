function newton_raphson_new()
    run('PR_Mansoori.m');
    
    initial_z_values = 0:10^-4:1;
    num_initial_values = length(initial_z_values);
    
    % Define pressure range
    P_range = 1:1:300;
    num_pressures = length(P_range);

    % Arrays to store data for all pressure points
    roots = zeros(num_initial_values, num_pressures);
    phi_values = zeros(num_initial_values, num_pressures);
    tau_values = zeros(num_initial_values, num_pressures);
    convergence_flags = zeros(num_initial_values, num_pressures);
    iteration_counts = zeros(num_initial_values, num_pressures);
    calculated_y1_values = zeros(num_initial_values, num_pressures);
    relative_differences = zeros(num_initial_values, num_pressures);

    T_fixed = 308;
    
    % Parameters for calculating sublimation pressure
    T_ref = 304.18;        % Reference temperature in K
    P_ref = 1e-6;          % Reference pressure in bar (assumed very low for pharmaceutical compounds)
    delta_H_vap = 56910;   % Enthalpy of vaporization in J/mol (56.91 kJ/mol)
    R_gas = 8.314;         % Gas constant in J/(mol·K)
    
    % Calculate sublimation pressure
    P_sublimation = P_ref * exp((delta_H_vap/R_gas)*(1/T_ref - 1/T_fixed));
    fprintf('Calculated sublimation pressure at %d K: %.4e bar\n', T_fixed, P_sublimation);
    
    % Calculate solid molar volume (MW/density)
    MW_drug = 598.5;       % Molecular weight of Ceftriaxone sodium in g/mol
    density_solid = 1.5;   % Assumed density in g/cm³ (typical for pharmaceutical solids)
    v_solid = MW_drug/density_solid; % Molar volume in cm³/mol
    fprintf('Assumed solid molar volume: %.2f cm³/mol\n', v_solid);
    
    % Convert to m³/mol for calculation
    v_solid_SI = v_solid * 1e-6; % Convert cm³/mol to m³/mol
    
    % Convert sublimation pressure to Pa for calculation
    P_sublimation_Pa = P_sublimation * 1e5;
    
    % Arrays to collect filtered results
    filtered_calculated_mole_fractions = [];
    filtered_initial_mole_fractions = [];
    filtered_rel_differences = [];
    filtered_pressures = [];
    filtered_temperatures = [];
    
    % Process each pressure value
    for p_idx = 1:num_pressures
        P_fixed = P_range(p_idx);
        
        % Skip very low pressures to avoid numerical issues
        if P_fixed < 5
            continue;
        end
        
        fprintf('\nProcessing pressure: %d bar (%d of %d)\n', P_fixed, p_idx, num_pressures);
        
        % Convert current pressure to Pa
        P_fixed_Pa = P_fixed * 1e5;

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
          
            tau_values(i, p_idx) = tau_1;

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
                roots(i, p_idx) = NaN;
                phi_values(i, p_idx) = NaN;
                calculated_y1_values(i, p_idx) = NaN;
                relative_differences(i, p_idx) = NaN;
                convergence_flags(i, p_idx) = 0;
                iteration_counts(i, p_idx) = max_iterations;
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
                    calculated_y1_values(i, p_idx) = NaN;  % Cannot calculate y_1 with phi = 0
                    relative_differences(i, p_idx) = NaN;
                elseif ln_phi > 700
                    phi = Inf;
                    calculated_y1_values(i, p_idx) = NaN;  % Cannot calculate y_1 with phi = Inf
                    relative_differences(i, p_idx) = NaN;
                else
                    phi = exp(ln_phi);
                    
                    % Now calculate y_1 using the formula:
                    % y_1 = (φ^Saturated_1 * P_1^Sublimation * exp[v_1^solid(P-P_1^Sublimation)/(RT)]) / (P * φ^supercritical_1)
                    % With φ^Saturated_1 = 1
                    
                    % Calculate Poynting factor
                    poynting_factor = exp((v_solid_SI * (P_fixed_Pa - P_sublimation_Pa)) / (R_gas * T_fixed));
                    
                    % Calculate new y_1
                    new_y_1 = (1 * P_sublimation * poynting_factor) / (P_fixed * phi);
                    calculated_y1_values(i, p_idx) = new_y_1;
                    
                    % Calculate relative difference
                    rel_diff = abs(initial_z_values(i) - new_y_1) / new_y_1;
                    relative_differences(i, p_idx) = rel_diff;
                    
                    % Add to filtered results if relative difference is < 5%
                    if rel_diff < 0.05
                        filtered_calculated_mole_fractions = [filtered_calculated_mole_fractions; new_y_1];
                        filtered_initial_mole_fractions = [filtered_initial_mole_fractions; y_1];
                        filtered_rel_differences = [filtered_rel_differences; rel_diff * 100]; % Convert to percentage
                        filtered_pressures = [filtered_pressures; P_fixed];
                        filtered_temperatures = [filtered_temperatures; T_fixed];
                    end
                end
                
                roots(i, p_idx) = root;
                phi_values(i, p_idx) = phi;
                convergence_flags(i, p_idx) = converged;
                iteration_counts(i, p_idx) = iterations;
            end
        end
        
        % Find and display the best match for this pressure
        valid_indices = ~isnan(calculated_y1_values(:, p_idx));
        if any(valid_indices)
            rel_diffs = relative_differences(:, p_idx);
            valid_rel_diffs = rel_diffs(valid_indices);
            valid_initial_y1 = initial_z_values(valid_indices);
            valid_calculated_y1 = calculated_y1_values(valid_indices, p_idx);
            
            [min_diff, min_idx] = min(valid_rel_diffs);
            
            fprintf('Best match at P = %d bar: Initial y_1 = %.6e, Calculated y_1 = %.6e, Relative Difference = %.6f%%\n', ...
                P_fixed, valid_initial_y1(min_idx), valid_calculated_y1(min_idx), min_diff*100);
        end
    end
    
    % Create the comparison table with filtered results
    comparison_table = table(...
        filtered_calculated_mole_fractions, ...
        filtered_initial_mole_fractions, ...
        filtered_rel_differences, ...
        filtered_pressures, ...
        filtered_temperatures, ...
        'VariableNames', {'CalculatedMoleFraction', 'InitialMoleFraction', 'RelativeDifferencePercent', 'Pressure', 'Temperature'});
    
    % Sort by pressure
    comparison_table = sortrows(comparison_table, 'Pressure');
    
    % Save to CSV
    writetable(comparison_table, 'mole_fraction_comparison.csv');
    fprintf('\nFiltered comparison data (rel diff < 5%%) saved to mole_fraction_comparison.csv\n');
    fprintf('Found %d points with relative difference < 5%%\n', height(comparison_table));
    
    % Create plot of best matching y_1 values vs pressure
    best_match_y1 = zeros(num_pressures, 1);
    best_match_calculated_y1 = zeros(num_pressures, 1);
    best_match_rel_diff = zeros(num_pressures, 1);
    
    for p_idx = 1:num_pressures
        valid_indices = ~isnan(calculated_y1_values(:, p_idx));
        if any(valid_indices)
            rel_diffs = relative_differences(:, p_idx);
            valid_rel_diffs = rel_diffs(valid_indices);
            valid_initial_y1 = initial_z_values(valid_indices);
            valid_calculated_y1 = calculated_y1_values(valid_indices, p_idx);
            
            [min_diff, min_idx] = min(valid_rel_diffs);
            
            best_match_y1(p_idx) = valid_initial_y1(min_idx);
            best_match_calculated_y1(p_idx) = valid_calculated_y1(min_idx);
            best_match_rel_diff(p_idx) = min_diff;
        else
            best_match_y1(p_idx) = NaN;
            best_match_calculated_y1(p_idx) = NaN;
            best_match_rel_diff(p_idx) = NaN;
        end
    end
    
    valid_p_indices = ~isnan(best_match_y1);
    
    figure;
    semilogy(P_range(valid_p_indices), best_match_calculated_y1(valid_p_indices), 'ro-', 'LineWidth', 2);
    hold on;
    semilogy(P_range(valid_p_indices), best_match_y1(valid_p_indices), 'b--', 'LineWidth', 1);
    grid on;
    xlabel('Pressure (bar)');
    ylabel('Mole Fraction (log scale)');
    title('Best Matching Initial y_1 vs. Calculated y_1 Across Pressure Range');
    legend('Calculated y_1', 'Initial y_1', 'Location', 'northwest');
    
    % Save full results to CSV (for reference)
    all_results = [];
    for p_idx = 1:num_pressures
        valid_indices = ~isnan(calculated_y1_values(:, p_idx));
        if any(valid_indices)
            P_current = P_range(p_idx);
            curr_results = [initial_z_values(valid_indices)', ...
                            roots(valid_indices, p_idx), ...
                            phi_values(valid_indices, p_idx), ... 
                            tau_values(valid_indices, p_idx), ...
                            calculated_y1_values(valid_indices, p_idx), ...
                            relative_differences(valid_indices, p_idx), ...
                            repmat(P_current, sum(valid_indices), 1), ...
                            repmat(T_fixed, sum(valid_indices), 1)];
            all_results = [all_results; curr_results];
        end
    end
    
    if ~isempty(all_results)
        results_table = array2table(all_results, 'VariableNames', ...
            {'InitialMoleFraction', 'Z', 'Phi', 'Tau', 'CalculatedMoleFraction', ...
             'RelativeDifference', 'Pressure', 'Temperature'});
        writetable(results_table, 'fugacity_results_all_pressures.csv');
        fprintf('All results saved to fugacity_results_all_pressures.csv\n');
    end
end
% '
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