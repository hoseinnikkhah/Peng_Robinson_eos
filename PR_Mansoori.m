% Genral info
Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;
T = [308 318 328 338]; % Temp range
P = [120 150 180 210 240 270]; % Pressure range
% The thermodynamic properties of Ceftriaxone sodium
T_c_drug = 304.18; % Critical Temp [K]
P_c_drug = 73.8; % Critical Pressure [bar]
omega_drug = 0.225;
% The thermodynamic properties of CO2
T_c_CO2 = 1149.5; % Critical Temp [K]
P_c_CO2 = 10.57; % Critical Pressure [bar]
omega_CO2 = 2.0964;
data_CO2 = zeros(3,24);
data_drug = zeros(3,24);
count = 1;
for i=1:4
    for j=1:6
        [a, b, d, A, B, D] = correlations(T_c_CO2,P_c_CO2,T(i),P(j),omega_CO2);
        data_CO2(:, count) = [A; B; D];
        [a, b, d, A, B, D] = correlations(T_c_CO2,P_c_drug,T(i),P(j),omega_drug);
        data_drug(:, count) = [A; B; D];
        count = count + 1;
    end
end
% coefficents
Z = zeros(1,24);
for i=1:24
    A_val = data_drug(1,i);
    B_val = data_drug(2,i);
    D_val = data_drug(3,i);
    
    c3 = 1;
    c2 = B_val - 1;
    c1 = D_val - 3*B_val^2 - 2*B_val + A_val - 2*sqrt(A_val*D_val);
    c0 = B_val^3 + B_val^2 - A_val*B_val - B_val*D_val + 2*B_val*sqrt(A_val*D_val);
    
    coeff = [c3, c2, c1, c0];
    root_values = roots(coeff);
    
    % Select only real roots
    real_roots = root_values(imag(root_values) == 0);
    
    % For gas phase calculations, typically the largest real root is used
    if ~isempty(real_roots)
        Z(1,i) = max(real_roots);
    else
        Z(1,i) = NaN; % Handle case where no real roots exist
    end
end