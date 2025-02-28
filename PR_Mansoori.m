% Genral info
Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;

T = [308 318 328 338];              % Temp range
P = [120 150 180 210 240 270];      % Pressure range
k_ij = [0.103 0.117 0.126 0.134];   % Binary interaction parameter for PR/KM EoS


% The thermodynamic properties of Ceftriaxone sodium
T_c_drug = 304.18;                  % Critical Temp [K]
P_c_drug = 73.8;                    % Critical Pressure [bar]
omega_drug = 0.225;

% The thermodynamic properties of CO2
T_c_CO2 = 1149.5;                   % Critical Temp [K]
P_c_CO2 = 10.57;                    % Critical Pressure [bar]
omega_CO2 = 2.0964;

data_CO2 = zeros(3,24);
data_drug = zeros(3,24);

data_CO2_cap = zeros(3,24);
data_drug_cap = zeros(3,24);

count = 1;

for i=1:4
    for j=1:6
        [a, b, d, A, B, D] = correlations(T_c_CO2,P_c_CO2,T(i),P(j),omega_CO2);
        data_CO2(:, count) = [a; b; d];
        data_CO2_cap(:, count) = [A; B; D];
        [a, b, d, A, B, D] = correlations(T_c_CO2,P_c_drug,T(i),P(j),omega_drug);
        data_drug(:, count) = [a; b; d];
        data_drug_cap(:, count) = [A; B; D];
        count = count + 1;
    end
end

data_CO2_3D = reshape(data_CO2, [6, 6, 4]);
data_drug_3D = reshape(data_drug, [6, 6, 4]);

data_CO2_3D_cap = reshape(data_CO2_cap, [6, 6, 4]);
data_drug_3D_cap = reshape(data_drug_cap, [6, 6, 4]);

a_ij = zeros(4,6);
for i=1:4
    for j=1:6
        factor = 6*(i-1) + j;
        a_ij(i,j) = sqrt(data_CO2(1,factor)*data_drug(1,factor))*k_ij(i);
    end
end

