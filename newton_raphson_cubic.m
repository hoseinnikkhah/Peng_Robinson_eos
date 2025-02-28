function [root, iterations, convergence] = newton_raphson_cubic(initial_guess, max_iter, tolerance, A, B, D)
    % Newton-Raphson method to solve the cubic equation:
    % Z^3 + (B−1)Z^2 + (D−3B^2−2B+A−2sqrt(AD))Z + (B^3 + B^2 − AB − BD + 2B*sqrt(AD)) = 0
    %
    % Inputs:
    % initial_guess - Starting value for Z
    % max_iter - Maximum number of iterations
    % tolerance - Convergence tolerance
    % A, B, D - Parameters for the equation
    %
    % Outputs:
    % root - Calculated root of the equation
    % iterations - Number of iterations performed
    % convergence - Boolean indicating whether solution converged

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
    
    fprintf('Starting Newton-Raphson method with initial guess z = %.6f\n', z);
    
    % Main iteration loop
    for i = 1:max_iter
        % Original function: f(z) = z^3 + coef_z2*z^2 + coef_z1*z + coef_z0
        f_z = z^3 + coef_z2*z^2 + coef_z1*z + coef_z0;
        
        % Derivative: f'(z) = 3z^2 + 2*coef_z2*z + coef_z1
        f_prime_z = 3*z^2 + 2*coef_z2*z + coef_z1;
        
        % Check if derivative is too close to zero
        if abs(f_prime_z) < 1e-10
            fprintf('Warning: Derivative close to zero, adjusting...\n');
            f_prime_z = sign(f_prime_z) * 1e-10;
            if f_prime_z == 0
                f_prime_z = 1e-10;  % Default to positive if sign is zero
            end
        end
        
        % Newton-Raphson update
        z_next = z - f_z / f_prime_z;
        
        % Calculate error for convergence check
        error = abs(z_next - z);
        
        fprintf('Iteration %d: z = %.10f, f(z) = %.10e, error = %.10e\n', i, z_next, f_z, error);
        
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
    if convergence
        fprintf('Solution converged after %d iterations: z = %.10f\n', iterations, root);
    else
        fprintf('Failed to converge after %d iterations. Last value: z = %.10f\n', iterations, root);
    end
end

% Example usage with parameter calculation
function [root] = solve_equation_with_calculated_parameters(initial_z, parameter_inputs)
    % This function first calculates parameters A, B, D based on provided inputs
    % Then solves the cubic equation using those parameters
    % 
    % Replace the parameter calculations below with your actual formulas
    % parameter_inputs should contain whatever values you need to calculate A, B, D
    
    % Example parameter calculations (REPLACE THESE with your actual formulas)
    A = calculate_A(parameter_inputs);
    B = calculate_B(parameter_inputs);
    D = calculate_D(parameter_inputs);
    
    fprintf('Calculated parameters: A = %.6f, B = %.6f, D = %.6f\n', A, B, D);
    
    % Solve using Newton-Raphson
    [root, iterations, convergence] = newton_raphson_cubic(initial_z, 100, 1e-10, A, B, D);
    
    % Check for convergence and perform additional operations if needed
    if ~convergence
        fprintf('Consider trying a different initial guess or increasing max iterations\n');
    end
end

% Example functions to calculate parameters (REPLACE THESE with your actual formulas)
function A = calculate_A(params)
    % Replace this with your formula for A
    A = params.some_value;  % Example placeholder
end

function B = calculate_B(params)
    % Replace this with your formula for B
    B = params.another_value;  % Example placeholder
end

function D = calculate_D(params)
    % Replace this with your formula for D
    D = params.third_value;  % Example placeholder
end