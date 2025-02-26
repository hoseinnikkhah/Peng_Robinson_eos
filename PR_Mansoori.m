% Genral info
Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;
T = [308 318 328 338];                  % Temp range
P = [120 150 180 210 240 270];          % Pressure range

% The thermodynamic properties of Ceftriaxone sodium
T_c_drug = 304.18;      % Critical Temp [K]
P_c_drug = 73.8;        % Critical Pressure [bar]
omega_drug = 0.225;

% The thermodynamic properties of CO2
T_c_CO2 = 1149.5;       % Critical Temp [K]
P_c_CO2 = 10.57;        % Critical Pressure [bar]
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
c3 = 1;
for i=1:24
    c2 = data_drug(2,i) - 1;
    c1 = data_drug(3,i) - 3*data_drug(2,i) + data_drug(1,i) - 2*sqrt(data_drug(1,i)*data_drug(3,i));
    c0 = data_drug(2,i)^3 + data_drug(2,i)^2 - data_drug(1,i)*data_drug(2,i) - data_drug(2,i)*data_drug(3,i) + 2*data_drug(2,i)*sqrt(data_drug(1,i)*data_drug(3,i));


