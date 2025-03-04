function [A, B, D] = dimensionless_factor(a_m, b_m, d_m, P, T)
    R = 8.314;    

    A = (a_m*P)/((R*T)^2);
    B = (b_m*P)/(R*T);
    D = (d_m*P)/(R*T);
end

