function [a, b, c ] = correlations_abd(T_c, P_c, T, m)
    R = 8.314;      [J/mol·K]

    a = 0.457235 * (((R^2)*(T_c^2))/P_c) * (1 + m*(1 - sqrt(T/T_c)))^2 * (1 + m)^2 
    
    