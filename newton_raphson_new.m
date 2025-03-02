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

        initial_guess = 0.2;
        max_iterations = 100;
        tolerance = 1e-10;
        
        [root, iterations, converged] = newton_raphson_cubic(initial_guess, max_iterations, tolerance, A, B, D);
        
        coeff1 = 2*(y_1*b_ij_3D(1,1,1) + y_2*b_ij_3D(1,2,1));
        coeff2 = (root + (1 + sqrt(2))*B)/(root + (1 - sqrt(2))*B);
        
        phi = exp((root-1)*((coeff1/b_m) - 1) - log(root - B) - (tau_1/(sqrt(2)*R*T_fixed*b_m))*log(coeff2));
        
        roots(i) = root;
        phi_values(i) = phi;
        convergence_flags(i) = converged;
        iteration_counts(i) = iterations;
        
        fprintf('y_1 = %.2f: Z = %.6f, phi = %.6e, Converged = %d, Iterations = %d\n', ...
                y_1, root, phi, converged, iterations);
    end
    
    figure;
    subplot(3,1,1);
    plot(initial_z_values, roots, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Z Value');
    title('Compressibility Factor (Z) vs. Mole Fraction');
    
    subplot(3,1,2);
    plot(initial_z_values, phi_values, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Phi Value');
    title('Fugacity Coefficient (Phi) vs. Mole Fraction');
    
    subplot(3,1,3);
    plot(initial_z_values, iteration_counts, 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Iterations');
    title('Convergence Iterations vs. Mole Fraction');
    
    non_converged = find(convergence_flags == 0);
    if ~isempty(non_converged)
        fprintf('\nWarning: Solution did not converge for the following mole fractions:\n');
        for i = 1:length(non_converged)
            fprintf('y_1 = %.2f\n', initial_z_values(non_converged(i)));
        end
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