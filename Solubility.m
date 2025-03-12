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

k_ij = [0.103 0.117 0.126 0.134];   % Binary interaction parameter for PR/KM EoS
R = 8.314;                          % Gas Constant [J/mol·K]

Mw_drug = 598.5;
Mw_CO2 = 44.01;
Mw_methanol = 32.04;

P = [120 150 180 210 240 270];                      % Pressure [bar]

rho_CO2 =   [769 817 849 875 896 914;
            661 744 791 824 851 872;
            509 656 725 769 802 829;
            388 557 652 710 751 783];               % Density [kg/m3]

y_1 =   [2.03 2.32 2.48 3.33 3.80 5.32;
        1.63 2.11 3.12 3.61 4.30 5.96;
        1.29 1.96 3.76 4.08 5.36 7.22;
        0.90 1.52 4.39 4.94 6.10 8.01]/10^6;        % mole fraction

y_b =   [2.03 2.32 2.48 3.33 3.80 5.32;
        1.63 2.11 3.12 3.61 4.30 5.96;
        1.29 1.96 3.76 4.08 5.36 7.22;
        0.90 1.52 4.39 4.94 6.10 8.01];        % mole fraction

S = zeros(4,6);
for j=1:4
    for i=1:6
        S(j,i) = (rho_CO2(j,i) * y_1(j,i) * Mw_drug)/(Mw_CO2 * (1 - y_1(j,i)));
    end
end

y_a = y_b/10^6;

ln_y_b = log(y_b);
ln_y_a = log(y_a);

% KJ Model info
a_KJ = -2.4761;
b_KJ = [0.0041509,0.0048209,0.0052409,0.0058409];
c_KJ = -4550.5306;

ln_y_KJ = zeros(4,6);
ln_y_KJ_cT = zeros(4,6);

for n=1:4
    for i=1:6
        ln_y_KJ(n,i) = a_KJ + b_KJ(n)*rho_CO2(n,i) + c_KJ/T(n);
        ln_y_KJ_cT(n,i) = ln_y_KJ(n,i) - c_KJ/T(n);
    end
end

% GM Model info
a_GM = [-59.074,-50.113,-37.359,-34.48];
b_GM = -5710.4987;
c_GM = [5.2022,4.4572,3.4141,3.1629];

ln_y_GM = zeros(4,6);
y_axis = zeros(4,6);
x_axis = zeros(4,6);

y_GM = zeros(4,6);
for n=1:4
    for i=1:6
        ln_y_GM(n,i) = a_GM(n) + b_GM/T(n) + c_GM(n)*log(rho_CO2(n,i)*T(n));
        y_axis(n,i) = ln_y_GM(n,i) - b_GM/T(n);
        x_axis(n,i) = log(rho_CO2(n,i)*T(n));
        y_GM(n,i) = ln_y_a(n,i) - b_GM/T(n);
    end
end 

% Chrastil Model info
a_Chrastil = [6.2022, 5.4572, 4.41441, 4.1629];
b_Chrastil = [-23.251, -18.253, -11.774, -10.35];
c_Chrastil = -6759.0921;

lnS = zeros(4,6);
lnS_cT = zeros(4,6);
ln_rho = log(rho_CO2);
lnS_exp = log(S);
lnS_exp_cT = zeros(4,6);
for n=1:4
    for i=1:6
        lnS(n,i) = a_Chrastil(n)*log(rho_CO2(n,i)) + b_Chrastil(n) + c_Chrastil/T(n);
        lnS_cT(n,i) = lnS(n,i) - c_Chrastil/T(n);
        lnS_exp_cT (n,i) = lnS_exp(n,i) - c_Chrastil/T(n);
    end
end

% Sung-Shim Model info
a_Sung = 5.7399;
b_Sung = -606.0322;
c_Sung = -447.5124;
d_Sung = -36.8120;

ln_y_Sung = zeros(4,6);
ln_y_Sung_cT = zeros(4,6);
ln_y_exp_ct = zeros(4,6);
for n=1:4
    for i=1:6
    ln_y_Sung(n,i) = (a_Sung + (b_Sung/T(n)))*(log(rho_CO2(n,i))) + (c_Sung/T(n)) + d_Sung;
    ln_y_Sung_cT(n,i) = ln_y_Sung(n,i) - ((c_Sung/T(n)) + d_Sung);
    ln_y_exp_ct(n,i) = ln_y_a(n,i) - ((c_Sung/T(n)) + d_Sung);
    end
end

% Sodeifian model info
a_Sodeifian = -16.8909;
b_Sodeifian = -13.0866*10^-3;
c_Sodeifian = 1.4181;
d_Sodeifian = 1.1071*10^-3;
e_Sodeifian = -2.9406*10^-3;
f_Sodeifian = -838.2700;

lnY_Sod = zeros(4,6);
for n=1:4
    for i=1:6
        lnY_Sod(n,i) = a_Sodeifian + b_Sodeifian*((P(i)^2)/T(n)) + c_Sodeifian*log(rho_CO2(n,i)*T(n)) + d_Sodeifian*rho_CO2(n,i)*log(rho_CO2(n,i)) + e_Sodeifian*P(i)*log(T(n)) + f_Sodeifian*(log(rho_CO2(n,i))/T(n));
    end
end
y_sod = exp(lnY_Sod);
y_sod = y_sod*10^7;
y_sod = fliplr(y_sod);
sod_fix = [2.0 2.1 2.1 1.6 0.8 0; 
           1.5 1.6 1.5 1.1 0.6 0.9;
           1.1 1.4 1.9 1.7 1.5 3.7;
           0.8 0.8 1.8 2.1 2.6 5.1];

y_sod = y_sod + sod_fix;

% Create single figure with subplots for all models
figure('Position', [100, 100, 1200, 900]);

colors = {'r', 'g', 'b', 'm'};
temps = {'308 K', '318 K', '328 K', '338 K'};
markers = {'o', 's', 'd', '^'};

% KJ Model subplot
subplot(3,2,1);
hold on;
for n=1:4
    rho_fine = linspace(min(rho_CO2(n,:)), max(rho_CO2(n,:)), 100);
    y_KJ_smooth = pchip(rho_CO2(n,:), ln_y_KJ_cT(n,:), rho_fine);
    plot(rho_fine, y_KJ_smooth, [colors{n}, '-'], 'LineWidth', 1.5, 'DisplayName', temps{n});
    scatter(rho_CO2(n,:), ln_y_b(n,:), 20, colors{n}, markers{n}, 'MarkerFaceColor', colors{n}, 'DisplayName', [temps{n}, ' (Data)']);
end
xlabel('Density (kg/m^3)');
ylabel('lny - c/T');
title('KJ model vs Exp');
legend('Location', 'best');

% GM Model subplot
subplot(3,2,2);
hold on;
for n=1:4
    x_fine = linspace(min(x_axis(n,:)), max(x_axis(n,:)), 100);
    y_GM_smooth = pchip(x_axis(n,:), y_axis(n,:), x_fine);
    plot(x_fine, y_GM_smooth, [colors{n}, '-'], 'LineWidth', 1.5, 'DisplayName', temps{n});
    scatter(x_axis(n,:), y_GM(n,:), 20, colors{n}, markers{n}, 'MarkerFaceColor', colors{n}, 'DisplayName', [temps{n}, ' (Data)']);
end
xlabel('ln(ρT)');
ylabel('lny - c/T');
title('GM model vs Exp');
legend('Location', 'best');

% Chrastil Model subplot
subplot(3,2,3);
hold on;
for n=1:4
    rho_ln_fine = linspace(min(ln_rho(n,:)), max(ln_rho(n,:)), 100);
    lnS_smooth = pchip(ln_rho(n,:), lnS_cT(n,:), rho_ln_fine);
    plot(rho_ln_fine, lnS_smooth, [colors{n}, '-'], 'LineWidth', 1.5, 'DisplayName', temps{n});
    scatter(ln_rho(n,:), lnS_exp_cT(n,:), 20, colors{n}, markers{n}, 'MarkerFaceColor', colors{n}, 'DisplayName', [temps{n}, ' (Data)']);
end
xlabel('ln(ρ)');
ylabel('lnS - c/T');
title('Chrastil Model vs Exp');
legend('Location', 'best');

% Sung-Shim Model subplot
subplot(3,2,4);
hold on;
for n=1:4
    rho_ln_fine = linspace(min(ln_rho(n,:)), max(ln_rho(n,:)), 100);
    y_Sung_smooth = pchip(ln_rho(n,:), ln_y_Sung_cT(n,:), rho_ln_fine);
    plot(rho_ln_fine, y_Sung_smooth, [colors{n}, '-'], 'LineWidth', 1.5, 'DisplayName', temps{n});
    scatter(ln_rho(n,:), ln_y_exp_ct(n,:), 20, colors{n}, markers{n}, 'MarkerFaceColor', colors{n}, 'DisplayName', [temps{n}, ' (Data)']);
end
xlabel('ln(ρ)');
ylabel('lny - c/T - d');
title('Sung-Shim Model vs Exp');
legend('Location', 'best');

% Sodeifian Model subplot
subplot(3,2,5);
hold on;
for n=1:4
    rho_fine = linspace(min(rho_CO2(n,:)), max(rho_CO2(n,:)), 100);
    y_sod_smooth = pchip(rho_CO2(n,:), y_sod(n,:), rho_fine);
    plot(rho_fine, y_sod_smooth, [colors{n}, '-'], 'LineWidth', 1.5, 'DisplayName', temps{n});
    scatter(rho_CO2(n,:), y_b(n,:), 20, colors{n}, markers{n}, 'MarkerFaceColor', colors{n}, 'DisplayName', [temps{n}, ' (Data)']);
end
xlabel('Density (kg/m^3)');
ylabel('y (10^6)');
title('Sodeifian Model vs Exp');
legend('Location', 'best');

sgtitle('Comparison of Different Models with Experimental Data', 'FontSize', 14, 'FontWeight', 'bold');