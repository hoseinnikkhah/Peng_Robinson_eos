function [a, b, d] = correlations_cap(T, P, [a; b; d])
    R = 8.314;      % [J/mol·K]

    A = (a*P)/((R*T)^2);
    B = (b*P)/(R*T);
    D = (d*P)/(R*T);
end

