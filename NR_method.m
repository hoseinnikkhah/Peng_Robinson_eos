% Main script to solve cubic equation with varying initial guesses
function solve_cubic_equation_varying_initial_guesses()
    % Run the Peng-Robinson EoS calculations
    run('PR_Mansoori.m');
    
    % Define the range of initial guesses (mole fractions)
    initial_z_values = 0:0.01:1;
    num_initial_values = length(initial_z_values);
    
    % Arrays to store results
    roots = zeros(size(initial_z_values));
    convergence_flags = zeros(size(initial_z_values));
    iteration_counts = zeros(size(initial_z_values));
    
    % Fixed temperature and pressure for testing
    T_fixed = 308;    % Assuming this matches T(1) in your code
    P_fixed = 120;    % Assuming this matches P(1) in your code
    
    % Process each initial guess
    for i = 1:num_initial_values
        y_1 = initial_z_values(i);
        y_2 = 1 - y_1;
        
        % Calculate mixture parameters based on mole fractions
        b_m = (y_1^2 * b_ij_3D(1,1,1)) + (2*y_1*y_2*b_ij_3D(1,2,1)) + (y_2^2 * b_ij_3D(2,2,1));
        d_m = (y_1^2 * d_ij_3D(1,1,1)) + (2*y_1*y_2*d_ij_3D(1,2,1)) + (y_2^2 * d_ij_3D(2,2,1));
        upper = (y_1^2 * a_ij_3D(1,1,1)^(2/3) * b_ij_3D(1,1,1)^(1/3)) + ...
                (2*y_1*y_2 * a_ij_3D(1,2,1)^(2/3) * b_ij_3D(1,2,1)^(1/3)) + ...
                (y_2^2 * a_ij_3D(2,2,1)^(2/3) * b_ij_3D(2,2,1)^(1/3));
        lower = (y_1^2 * b_ij_3D(1,1,1)) + (2*y_1*y_2*b_ij_3D(1,2,1)) + (y_2^2 * b_ij_3D(2,2,1));
        a_m = (upper^1.5) / sqrt(lower);
        
        % Calculate A, B, D parameters
        [A, B, D] = correlations_cap(T_fixed, P_fixed, a_m, b_m, d_m);
        
        % Use Newton-Raphson to solve the cubic equation
        initial_guess = 0.5;  % You can adjust this starting value if needed
        max_iterations = 100;
        tolerance = 1e-10;
        
        [root, iterations, converged] = newton_raphson_cubic(initial_guess, max_iterations, tolerance, A, B, D);
        
        % Store results
        roots(i) = root;
        convergence_flags(i) = converged;
        iteration_counts(i) = iterations;
        
        % Display progress
        fprintf('Initial mole fraction y_1 = %.2f: Root = %.10f, Converged = %d, Iterations = %d\n', ...
                y_1, root, converged, iterations);
    end
    
    % Plot results
    figure;
    subplot(2,1,1);
    plot(initial_z_values, roots, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Root Value (Z)');
    title('Root Values vs. Mole Fraction');
    
    subplot(2,1,2);
    plot(initial_z_values, iteration_counts, 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Mole Fraction (y_1)');
    ylabel('Iterations');
    title('Convergence Iterations vs. Mole Fraction');
    
    % Identify non-converged solutions if any
    non_converged = find(convergence_flags == 0);
    if ~isempty(non_converged)
        fprintf('\nWarning: Solution did not converge for the following mole fractions:\n');
        for i = 1:length(non_converged)
            fprintf('y_1 = %.2f\n', initial_z_values(non_converged(i)));
        end
    else
        fprintf('\nAll solutions converged successfully.\n');
    end
    
    % Save results to a file
    results_table = table(initial_z_values', roots', convergence_flags', iteration_counts', ...
                          'VariableNames', {'MoleFraction', 'Root', 'Converged', 'Iterations'});
    writetable(results_table, 'cubic_equation_roots.csv');
    fprintf('Results saved to cubic_equation_roots.csv\n');
end

function [root, iterations, convergence] = newton_raphson_cubic(initial_guess, max_iter, tolerance, A, B, D)
    % Newton-Raphson method to solve the cubic equation:
    % Z^3 + (B−1)Z^2 + (D−3B^2−2B+A−2sqrt(AD))Z + (B^3 + B^2 − AB − BD + 2B*sqrt(AD)) = 0
    
    % Initialize variables
    z = initial_guess;
    iterations = 0;
    convergence = false;
    
    % Calculate sqrt(AD) once to avoid repeated calculation
    sqrt_AD = sqrt(A*D);
    
    % Coefficient of Z^2
    coef_z2 = B - 1;
    
    % Coefficient of Z^1
    coef_z1 = D - 3*B^2 - 2*B + A - 2*sqrt_AD;
    
    % Constant term
    coef_z0 = B^3 + B^2 - A*B - B*D + 2*B*sqrt_AD;
    
    % Main iteration loop
    for i = 1:max_iter
        % Original function: f(z) = z^3 + coef_z2*z^2 + coef_z1*z + coef_z0
        f_z = z^3 + coef_z2*z^2 + coef_z1*z + coef_z0;
        
        % Derivative: f'(z) = 3z^2 + 2*coef_z2*z + coef_z1
        f_prime_z = 3*z^2 + 2*coef_z2*z + coef_z1;
        
        % Check if derivative is too close to zero
        if abs(f_prime_z) < 1e-10
            % If derivative is near zero, make a small adjustment to avoid division by zero
            f_prime_z = sign(f_prime_z) * 1e-10;
            if f_prime_z == 0
                f_prime_z = 1e-10;  % Default to positive if sign is zero
            end
        end
        
        % Newton-Raphson update
        z_next = z - f_z / f_prime_z;
        
        % Calculate error for convergence check
        error = abs(z_next - z);
        
        % Update z for next iteration
        z = z_next;
        iterations = i;
        
        % Check for convergence
        if error < tolerance
            convergence = true;
            break;
        end
    end
    
    % Final output
    root = z;
end