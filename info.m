T_c = 253;
Mw_drug = 598.5;
Mw_CO2 = 44.01;

P = [120 150 180 210 240 270];                      % Pressure [bar]
rho_CO2 = [769 817 849 875 896 914];                % Density [kg/m3]
y_1 = [2.03 2.32 2.48 3.33 3.80 5.32]/10^6;         % mole fraction

S = linspace(0,1,6);
for i=1:6
    S(i) = (rho_CO2(i) * y_1(i) * Mw_drug)/(Mw_CO2 * (1 - y_1(i)));

end

% C_drug = [0.021 0.026 0.029 0.039 0.046 0.066];     % solubility or concentration [g/L]

