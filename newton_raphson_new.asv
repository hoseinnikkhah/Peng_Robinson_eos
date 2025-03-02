function newton_raphson_new()
    run('PR_Mansoori.m');
    
    initial_z_values = 0:10^-4:1;
    num_initial_values = length(initial_z_values);

    roots = zeros(size(initial_z_values));
    phi_values = zeros(size(initial_z_values));
    tau_values = zeros(size(initial_z_values));
    convergence_flags = zeros(size(initial_z_values));
    iteration_counts = zeros(size(initial_z_values));

    T_fixed = 308;
    P_fixed = 120;

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
            fprintf('Warning: Could not find any valid root for y_1 = %.4f\n', y_1);
            roots(i) = NaN;
            phi_values(i) = NaN;
            convergence_flags(i) = 0;
            iteration_counts(i) = max_iterations;
        else
            root = best_root;
            iterations = best_iterations;
            converged = best_converged;
            
            % If root is less than or equal to B, adjust it slightly
            if root <= B
                root = B + 1e-6;  % Ensure we are slightly above B for log(Z-B)
                fprintf('Warning: Using adjusted root (Z = %.6f) for y_1 = %.4f\n', root, y_1);
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
                fprintf('Warning: ln(phi) too negative (%.2e) for y_1 = %.4f, setting phi = 0\n', ln_phi, y_1);
            elseif ln_phi > 700
                phi = Inf;
                fprintf('Warning: ln(phi) too positive (%.2e) for y_1 = %.4f, setting phi = Inf\n', ln_phi, y_1);
            else
                phi = exp(ln_phi);
            end
            
            roots(i) = root;
            phi_values(i) = phi;
            convergence_flags(i) = converged;
            iteration_counts(i) = iterations;
            
            fprintf('y_1 = %.4f: Z = %.6f, phi = %.6e, ln(phi) = %.6e, Converged = %d, Iterations = %d\n', ...
                    y_1, root, phi, ln_phi, converged, iterations);
            
            % Add debugging for problematic cases
            if phi == 0 || phi == Inf || isnan(phi)
                fprintf('  Debug: term1 = %.6e, term2 = %.6e, term3 = %.6e\n', term1, term2, term3);
                fprintf('  Debug: tau_1 = %.6e, coeff1/b_m = %.6e, coeff2 = %.6e\n', tau_1, coeff1/b_m, coeff2);
            end
        end
    end
    
    % Filter out any NaN values for plotting
    valid_indices = ~isnan(roots);
    
    figure;
    subplot(3,1,1);
    plot(initial_z_values(valid_indices), roots(valid_indices), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Z Value');
    title('Compressibility Factor (Z) vs. Mole Fraction');
    
    subplot(3,1,2);
    semilogy(initial_z_values(valid_indices), phi_values(valid_indices), 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Phi Value (log scale)');
    title('Fugacity Coefficient (Phi) vs. Mole Fraction');
    
    subplot(3,1,3);
    plot(initial_z_values(valid_indices), tau_values(valid_indices), 'm-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Tau Value');
    title('Tau vs. Mole Fraction');
    
    non_converged = find(convergence_flags == 0);
    if ~isempty(non_converged)
        fprintf('\nWarning: Solution did not converge for %d mole fractions\n', length(non_converged));
    else
        fprintf('\nAll solutions converged successfully.\n');
    end
    
    results_table = table(initial_z_values', roots', phi_values', tau_values', convergence_flags', iteration_counts', ...
                          'VariableNames', {'MoleFraction', 'Z', 'Phi', 'Tau', 'Converged', 'Iterations'});
    writetable(results_table, 'fugacity_results.csv');
    fprintf('Results saved to fugacity_results.csv\n');
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