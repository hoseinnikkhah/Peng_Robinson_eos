% SC CO2 Thermo info
T_c_CO2 = 304.18;       % Critical Temp [Kelvin]
P_c_CO2 = 73.8;         % Critical Pressure [bar]
omega_CO2 = 0.225;      % Omega [-]

% Drug Thermo info
T_c_drug = 1149.5;      % Critical Temp [Kelvin]
P_c_drug = 10.57;       % Critical Pressure [bar]
omega_drug = 2.0964;    % Omega [-]

% Genral info
T = [308 318 328 338];              % Temp range
P = 1:1:300;                        % Pressure range from 1 to 300 bar with step of 1
k_ij = [0.103 0.117 0.126 0.134];   % Binary interaction parameter for PR/KM EoS
R = 8.314;                          % Gas Constant [J/mol·K]

Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;

correlations_CO2 = zeros(4,3);
correlations_drug = zeros(4,3);

for i=1:4
    [a, b, d] = correlations(T_c_CO2, P_c_CO2, T(i), omega_CO2);
    correlations_CO2(i,:) = [a,b,d];
    [a, b, d] = correlations(T_c_drug, P_c_drug, T(i), omega_drug);
    correlations_drug(i,:) = [a,b,d];
end


% a_ij calculations
% drug = 1 and CO2 = 2

a_ij = zeros(2,2,4);
b_ij = zeros(2,2,4);
d_ij = zeros(2,2,4);
for temp=1:4
    a_ij(1,1,temp) = correlations_drug(temp,1);
    a_ij(2,2,temp) = correlations_CO2(temp,1);
    a_ij(1,2,temp) = sqrt(correlations_drug(temp,1)*correlations_CO2(temp,1))*(1 - k_ij(1,temp));
    a_ij(2,1,temp) = sqrt(correlations_CO2(temp,1)*correlations_drug(temp,1))*(1 - k_ij(1,temp));

    b_ij(1,1,temp) = correlations_drug(temp,2);
    b_ij(2,2,temp) = correlations_CO2(temp,2);
    b_ij(1,2,temp) = (correlations_drug(temp,2) + correlations_CO2(temp,2))/2;
    b_ij(2,1,temp) = (correlations_CO2(temp,2) + correlations_drug(temp,2))/2;

    d_ij(1,1,temp) = correlations_drug(temp,3);
    d_ij(2,2,temp) = correlations_CO2(temp,3);
    d_ij(1,2,temp) = (correlations_drug(temp,3) + correlations_CO2(temp,3))/2;
    d_ij(2,1,temp) = (correlations_CO2(temp,3) + correlations_drug(temp,3))/2;
end

y_1 = 0:10^-8:0.999*10^-5;      % Drug mole fraction
y_2 = 1 - y_1;                  % CO2 mole fraction
mole_length = length(y_1);

mixed = zeros(3,mole_length,4);
for n=1:4
    for i=1:mole_length
        [a_mm, b_mm, d_mm] = mixing_rules(y_1(i), y_2(i), a_ij, b_ij, d_ij, n);
        mixed(:,i,n) = [a_mm; b_mm; d_mm];
    end
end

count = 1;
ABD_308 = zeros(3,mole_length,300);
ABD_318 = zeros(3,mole_length,300);
ABD_328 = zeros(3,mole_length,300);
ABD_338 = zeros(3,mole_length,300);

for n=1:4
    for p=1:300
        for i=1:mole_length
            [A, B, D] = dimensionless_factor(mixed(1,i,n), mixed(2,i,n), mixed(3,i,n), P(p), T(n));
            if n == 1
                ABD_308(:,i,p) = [A; B; D];
            elseif n == 2
                ABD_318(:,i,p) = [A; B; D];
            elseif n == 3
                ABD_328(:,i,p) = [A; B; D];
            elseif n == 4
                ABD_338(:,i,p) = [A; B; D];                                
            end
        end
    end
end



ABD_all = zeros(3,mole_length,300,4);
ABD_all(:,:,:,1) = ABD_308;
ABD_all(:,:,:,2) = ABD_318;
ABD_all(:,:,:,3) = ABD_328;
ABD_all(:,:,:,4) = ABD_338;

Z_roots_all = zeros(300,mole_length,4);

for temp=1:4
    for p=1:300
        for i=1:mole_length
            A = ABD_all(1,i,p,temp);
            B = ABD_all(2,i,p,temp);
            D = ABD_all(3,i,p,temp);

            C1 = 1;
            C2 = B-1;
            C3 = (D-3*B^2-2*B+A-2*sqrt(A*D));
            C4 = (B^3+B^2-A*B-B*D+2*B*sqrt(A*D));

            coefficients = [C1, C2, C3, C4];
            all_roots = roots(coefficients);

            real_roots = all_roots(abs(imag(all_roots)) < 1e-10);
            real_roots = real(real_roots);

            if length(real_roots) == 1
                Z = real_roots;
            elseif length(real_roots) == 3
                valid_roots = real_roots(real_roots > B);
    
                if isempty(valid_roots)
                    Z = max(real_roots);
                elseif length(valid_roots) == 1
                    Z = valid_roots;
                else
                    Z = max(valid_roots);
                end
            else
                Z = max(real_roots(real_roots > B));
            end

            Z_roots_all(p,i,temp) = Z;
        end
    end
end

tau_all = zeros(300,mole_length,4);
% This tau calculations and storing process
tau_308 = zeros(300,mole_length);
tau_318 = zeros(300,mole_length);
tau_328 = zeros(300,mole_length);
tau_338 = zeros(300,mole_length);

for temp = 1:4
    for p=1:300
        for i=1:mole_length
            [tau] = tau_calculation(y_1(i), y_2(i), mixed(1,i,temp), mixed(2,i,temp), mixed(3,i,temp), a_ij(1,1,temp), a_ij(1,2,temp), b_ij(1,1,temp), b_ij(1,2,temp), d_ij(1,1,temp), d_ij(1,2,temp), temp);
            if temp == 1
                tau_308(p,i) = tau;
            elseif temp == 2
                tau_318(p,i) = tau;
            elseif temp == 3
                tau_328(p,i) = tau;
            elseif temp == 4
                tau_338(p,i) = tau;
            end
        end
    end
end

tau_all(:,:,1) = tau_308;
tau_all(:,:,2) = tau_318;
tau_all(:,:,3) = tau_328;
tau_all(:,:,4) = tau_338;

phi = zeros(300,mole_length,4);

for temp=1:4
    for p=1:300
        for i=1:mole_length
            answer = phi_calculation(y_1(i), y_2(i), mixed(2,i,temp), b_ij(1,1,temp), b_ij(1,2,temp), T(temp), Z_roots_all(p,i,temp), ABD_all(2,i,p,n), tau_all(p,i,temp));
            phi(p,i,temp) = answer;
        end
    end
end

phi = log(phi);

% PSub_1 calculations
PSub_1 = zeros(1,4);
deltaH_vap = 56.91;             % [kJ/mol]
DeltaH_vap = deltaH_vap * 10^3; % [J/mol]
T_ref = 304.18;                 % [K]
P_ref = 10^-6;                  % [bar]
for temp=1:4
    PSub_1(temp) = P_ref * exp((DeltaH_vap/R)*(1/T_ref - 1/T(temp)));
end

phi_sat = 1;

rho_drug = 1.5; % [g/cm^3]
vSolid_1 = Mw_drug/rho_drug;
vSolid_1 = vSolid_1 * 10^-6; % [m^3/mol]

new_y_1 = zeros(300,mole_length,4);

for temp=1:4
    for p=1:300
        for i=1:mole_length
            EXP = (vSolid_1*(p-PSub_1(temp)))/(R*T(temp));
            new_y_1(p,i,temp) = ((phi_sat*PSub_1(temp))*exp(EXP))/(p*phi(p,i,temp));
        end
    end
end

difference_matrix = zeros(300,mole_length,4);

for temp=1:4
    for p=1:300
        for i=1:mole_length
            true_value = difference(new_y_1(p,i,temp), y_1(i), 5);
            difference_matrix(p,i,temp) = true_value;
        end
    end
end


extracted_values = cell(300, 4);  % [pressure, temperature]
extracted_indices = cell(300, 4); 


for temp = 1:4
    for p = 1:300
        
        valid_indices = find(~isnan(difference_matrix(p, :, temp)));
        valid_values = difference_matrix(p, valid_indices, temp);
        
        extracted_values{p, temp} = valid_values;
        extracted_indices{p, temp} = valid_indices;
    end
end


p_example = 238;
temp_example = 1;
valid_mole_fractions_at_p238_t1 = extracted_values{p_example, temp_example};

fprintf('Number of valid values at each pressure for temperature 1:\n');
for p = 1:300
    num_values = length(extracted_values{p, 1});
    if num_values > 0
        fprintf('Pressure %d: %d valid values\n', p, num_values);
    end
end