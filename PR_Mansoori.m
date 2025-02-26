% Genral info
Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;
T = 338;
P = 270;

% The thermodynamic properties of Ceftriaxone sodium
T_c_drug = 304.18;      % Temp [K]
P_c_drug = 73.8;        % Pressure [bar]
omega_drug = 0.225;

% The thermodynamic properties of CO2
T_c_CO2 = 1149.5;       % Temp [K]
P_c_CO2 = 10.57;        % Pressure [bar]
omega_CO2 = 2.0964;

[a, b, d, A, B, D] = correlations(T_c_CO2,P_c_CO2,T,P,omega_CO2);