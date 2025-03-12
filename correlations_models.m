clear;
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

rho_CO2 =   [769 817 849 875 896 914;
            661 744 791 824 851 872;
            509 656 725 769 802 829;
            388 557 652 710 751 783];               % Density [kg/m3]

y_b =   [2.03 2.32 2.48 3.33 3.80 5.32;
        1.63 2.11 3.12 3.61 4.30 5.96;
        1.29 1.96 3.76 4.08 5.36 7.22;
        0.90 1.52 4.39 4.94 6.10 8.01];        % mole fraction


% KJ Model info
a_KJ = -2.4761;
b_KJ = 5.2409*10^-3;
c_KJ = -4550.5306;

ln_y_KJ = zeros(4,6);
ln_y_KJ_cT = zeros(4,6);
ln_y_b = log(y_b);
for n=1:4
    for i=1:6
        ln_y_KJ(n,i) = a_KJ + b_KJ*rho_CO2(n,i) + c_KJ/T(n);
        ln_y_KJ_cT(n,i) = ln_y_KJ(n,i) - c_KJ/T(n);

    end
end

figure(1);
hold on;

% Plot lines
plot(rho_CO2(1,:), ln_y_KJ_cT(1,:), 'DisplayName', '308 K');
plot(rho_CO2(2,:), ln_y_KJ_cT(2,:), 'DisplayName', '318 K');
plot(rho_CO2(3,:), ln_y_KJ_cT(3,:), 'DisplayName', '328 K');
plot(rho_CO2(4,:), ln_y_KJ_cT(4,:), 'DisplayName', '338 K');

% Scatter points
scatter(rho_CO2(1,:), ln_y_b(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r', 'DisplayName', '308 K (Data)');
scatter(rho_CO2(2,:), ln_y_b(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g', 'DisplayName', '318 K (Data)');
scatter(rho_CO2(3,:), ln_y_b(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b', 'DisplayName', '328 K (Data)');
scatter(rho_CO2(4,:), ln_y_b(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm', 'DisplayName', '338 K (Data)');

% Axis labels
xlabel('Density (kg/m^3)');
ylabel('lny - c/T');

% Add legend
legend('Location', 'best');


% GM Model info
a_GM = 1.8309;
b_GM = -5710.4987;
c_GM = 3.8342;

% Chrastil Model info
a_Chrastil = 4.5889;
b_Chrastil = -12.7501;
c_Chrastil = -6759.0921;

% Bartle et al. Model info
a_Bartle = 13.2886;
b_Bartle = 9.0827;
c_Bartle = -6845.6966;

% Sung-Shim Model info
a_Sung = 5.7399;
b_Sung = -606.0322;
c_Sung = -447.5124;
d_Sung = -36.8120;

% Bian et al model info
a_Bian = 4.2910;
b_Bian = -4806.5989;
c_Bian = 0.0355;
d_Bian = -1.0456;
e_Bian = -9.5861*10^4;

% Sodeifian model info
a_Sodeifian = -16.8909;
b_Sodeifian = -13.0866*10^3;
c_Sodeifian = 1.4181;
d_Sodeifian = 1.1071*10^3;
e_Sodeifian = -2.9406*10^3;
f_Sodeifian = -838.2700;

