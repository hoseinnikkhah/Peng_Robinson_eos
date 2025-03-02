% Genral info
Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;

T = [308 318 328 338];              % Temp range
P = 1:1:300;                        % Pressure range from 1 to 300 bar with step of 1

k_ij = [0.103 0.117 0.126 0.134];   % Binary interaction parameter for PR/KM EoS
R = 8.314;

% The thermodynamic properties of Ceftriaxone sodium
T_c_drug = 304.18;                  % Critical Temp [K]
P_c_drug = 73.8;                    % Critical Pressure [bar]
omega_drug = 0.225;

% The thermodynamic properties of CO2
T_c_CO2 = 1149.5;                   % Critical Temp [K]
P_c_CO2 = 10.57;                    % Critical Pressure [bar]
omega_CO2 = 2.0964;

% Number of pressure points
num_pressures = length(P);

data_CO2 = zeros(3,4);
data_drug = zeros(3,4);

data_CO2_cap = zeros(3, 4*num_pressures);  % Now 3 × (4*300) = 3 × 1200
data_drug_cap = zeros(3, 4*num_pressures); % Now 3 × (4*300) = 3 × 1200

count = 1;

for i=1:4
    [a, b, d] = correlations(T_c_CO2, P_c_CO2, T(i), omega_CO2);
    data_CO2(:, i) = [a; b; d];                     % Checked and it works
    [a, b, d] = correlations(T_c_drug, P_c_drug, T(i), omega_drug);
    data_drug(:, i) = [a; b; d];                    % Checked and it works
end

for i=1:4
    for j=1:num_pressures
        [A, B, D] = correlations_cap(T(i), P(j), data_CO2(1,i), data_CO2(2,i), data_CO2(3,i));
        data_CO2_cap(:, count) = [A; B; D];         % Checked and it works
        [A, B, D] = correlations_cap(T(i), P(j), data_drug(1,i), data_drug(2,i), data_drug(3,i));
        data_drug_cap(:, count) = [A; B; D];        % Checked and it works
        count = count + 1;
    end
end

a_3D = zeros(2, 1, 4);
b_3D = zeros(2, 1, 4);
d_3D = zeros(2, 1, 4);

data_CO2_3D = reshape(data_CO2_cap, [3, num_pressures, 4]);  % Now [3, 300, 4]
data_drug_3D = reshape(data_drug_cap, [3, num_pressures, 4]); % Now [3, 300, 4]

for i = 1:4                                         % Checked and it works
    a_3D(1, 1, i) = data_CO2(1, i); 
    a_3D(2, 1, i) = data_drug(1, i);  
    b_3D(1, 1, i) = data_CO2(2, i); 
    b_3D(2, 1, i) = data_drug(2, i); 
    d_3D(1, 1, i) = data_CO2(3, i); 
    d_3D(2, 1, i) = data_drug(3, i); 
end

a_ij_3D = zeros(2, 2, 4);
b_ij_3D = zeros(2, 2, 4);
d_ij_3D = zeros(2, 2, 4);

for temp=1:4
    for i=1:2
        for j=1:2
            if i == j
                a_ij_3D(i,j,temp) = a_3D(i,1,temp);
                b_ij_3D(i,j,temp) = b_3D(i,1,temp);
                d_ij_3D(i,j,temp) = d_3D(i,1,temp);
            else
                a_ij_3D(i,j,temp) = sqrt(a_3D(i,1,temp)*a_3D(j,1,temp))*(1-k_ij(1,temp));
                b_ij_3D(i,j,temp) = (b_3D(i,1,temp) + b_3D(j,1,temp))/2;
                d_ij_3D(i,j,temp) = (d_3D(i,1,temp) + d_3D(j,1,temp))/2;
            end
        end
    end
end