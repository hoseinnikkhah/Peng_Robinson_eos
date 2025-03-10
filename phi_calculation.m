function [phi] = phi_calculation(y_1, y_2, b_m, b_11, b_12, T, Z, B, tau)
    R = 8.314;

    % Ensure Z > B with a safe margin
    if Z <= B || Z - B < 1e-8
        Z = B + 1e-6;
    end
    
    term1 = (Z-1)*((2*((y_1*b_11) + (y_2*b_12)))/b_m - 1);
    
    % Protect against ln(0) or ln(negative)
    if Z - B <= 0
        term2 = -log(1e-10);  % Extremely large negative number
    else
        term2 = -log(Z - B);
    end
    
    term3_numerator = Z + (1 + sqrt(2))*B;
    term3_denominator = Z + (1 - sqrt(2))*B;
    
    % Protect against division by zero
    if term3_denominator <= 0
        term3_denominator = 1e-10;
    end
    
    term3 = -(tau/(sqrt(2)*R*T*b_m))*log(term3_numerator/term3_denominator);
    
    ln_phi = term1 + term2 + term3;
    
    % Protect against overflow in exp()
    if ln_phi > 700
        phi = 1e300;  % Large but finite number
    elseif ln_phi < -700
        phi = 1e-300; % Small but non-zero number
    else
        phi = exp(ln_phi);
    end
end