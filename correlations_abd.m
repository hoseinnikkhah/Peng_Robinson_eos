function [a, b, c ] = correlations_abd(T_c, P_c, T, omega)
    R = 8.314;      [J/mol·K]

    m = 0.37464 + 1.5422*omega - 0.26992*(omega^2);

    a = 0.457235 * (((R^2)*(T_c^2))/P_c) * (1 + m*(1 - sqrt(T/T_c)))^2 * (1 + m)^2;
    b = 0.077796 * ((R*T_c)/P_c); 
    d = 0.457235 * (((R^2)*(T_c^2))/P_c) * (1 + m*(1 - sqrt(T/T_c)))^2 * ((m^2)/(R*T_c));
    
    