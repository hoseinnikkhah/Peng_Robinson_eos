function [tau] = tau_calculation(y_1, y_2, a_m, b_m, d_m, a_11, a_12, b_11, b_12, d_11, d_12, n)
    R = 8.314;
    if n == 1
        T = 308;
    elseif n == 2
        T = 318;
    elseif n == 3
        T = 328;
    elseif n == 4
        T = 338;
    end
    term1 = a_m + R*T*d_m;
    term2 = sqrt(a_m*b_m*R*T);
    term3 = ((y_1*b_11) + (y_2*b_12))/b_m;
    term4 = ((y_1*a_11) + (y_2*a_12));
    term5 = sqrt((R*T*d_m)/a_m);
    term6 = ((y_1*d_11) + (y_2*d_12));
    term7 = sqrt((R*T*a_m)/d_m);
    tau = term1 - 2*term2*(0.5 - term3) + term4*(1 - term5) + term6*((R*T) - term7);
end

