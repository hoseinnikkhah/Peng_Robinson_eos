function [A, B, D] = dimensionless_factor(mixing_correlations, P, T, n, i)
    R = 8.314;      % [J/mol·K]

    A = (mixing_correlations(1,i,n)*P)/((R*T)^2);
    B = (mixing_correlations(2,i,n)*P)/(R*T);
    D = (mixing_correlations(3,i,n)*P)/(R*T);
end

