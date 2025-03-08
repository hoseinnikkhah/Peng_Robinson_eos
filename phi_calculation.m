function [phi] = phi_calculation(y_1, y_2, b_m, T, Z, B)
    R = 8.314;

    term1 = ((2*((y_1*b_11) + (y_2*b_12)))/b_m) - 1;
    term2 = tau/(sqrt(2)*R*T*b_m);
    term3 = Z + (1 + sqrt(2))*B;
    term4 = Z + (1 - sqrt(2))*B;

    ln_phi = (Z-1)*term1 - log(Z - B) - term2*log(term3/term4);
    phi = exp(ln_phi);
end