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
% S = S*100;
% C_drug = [0.021 0.026 0.029 0.039 0.046 0.066];     % solubility or concentration [g/L]


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

xlabel('Density (kg/m^3)');
ylabel('lny - c/T');

legend('Location', 'best');
title('KJ model vs Exp');



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

figure(2);
hold on;

% Plot lines
plot(x_axis(1,:), y_axis(1,:), 'DisplayName', '308 K');
plot(x_axis(2,:), y_axis(2,:), 'DisplayName', '318 K');
plot(x_axis(3,:), y_axis(3,:), 'DisplayName', '328 K');
plot(x_axis(4,:), y_axis(4,:), 'DisplayName', '338 K');

% Scatter points
scatter(x_axis(1,:), y_GM(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r', 'DisplayName', '308 K (Data)');
scatter(x_axis(2,:), y_GM(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g', 'DisplayName', '318 K (Data)');
scatter(x_axis(3,:), y_GM(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b', 'DisplayName', '328 K (Data)');
scatter(x_axis(4,:), y_GM(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm', 'DisplayName', '338 K (Data)');

xlabel('Density (kg/m^3)');
ylabel('lny - c/T');

legend('Location', 'best');
title('GM model vs Exp');


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


figure(3);
hold on;

% Plot lines
plot(ln_rho(1,:), lnS_cT(1,:), 'DisplayName', '308 K');
plot(ln_rho(2,:), lnS_cT(2,:), 'DisplayName', '318 K');
plot(ln_rho(3,:), lnS_cT(3,:), 'DisplayName', '328 K');
plot(ln_rho(4,:), lnS_cT(4,:), 'DisplayName', '338 K');

% Scatter points
scatter(ln_rho(1,:), lnS_exp_cT(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r', 'DisplayName', '308 K (Data)');
scatter(ln_rho(2,:), lnS_exp_cT(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g', 'DisplayName', '318 K (Data)');
scatter(ln_rho(3,:), lnS_exp_cT(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b', 'DisplayName', '328 K (Data)');
scatter(ln_rho(4,:), lnS_exp_cT(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm', 'DisplayName', '338 K (Data)');

xlabel('Density (kg/m^3)');
ylabel('lny - c/T');

legend('Location', 'best');
title('Chrastil Model vs Exp');

% Bartle et al. Model info
% This model has an issue in the paper
a_Bartle = 13.2886;
b_Bartle = 9.0827;
c_Bartle = -6845.6966;


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

figure(4);
hold on;

% Plot lines
plot(ln_rho(1,:), ln_y_Sung_cT(1,:), 'DisplayName', '308 K');
plot(ln_rho(2,:), ln_y_Sung_cT(2,:), 'DisplayName', '318 K');
plot(ln_rho(3,:), ln_y_Sung_cT(3,:), 'DisplayName', '328 K');
plot(ln_rho(4,:), ln_y_Sung_cT(4,:), 'DisplayName', '338 K');

% Scatter points
scatter(ln_rho(1,:), ln_y_exp_ct(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r', 'DisplayName', '308 K (Data)');
scatter(ln_rho(2,:), ln_y_exp_ct(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g', 'DisplayName', '318 K (Data)');
scatter(ln_rho(3,:), ln_y_exp_ct(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b', 'DisplayName', '328 K (Data)');
scatter(ln_rho(4,:), ln_y_exp_ct(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm', 'DisplayName', '338 K (Data)');

xlabel('Density (kg/m^3)');
ylabel('lny - c/T');

legend('Location', 'best');
title('Sung-Shim Model vs Exp');



% Bian et al model info
a_Bian = 4.2910;
b_Bian = -4806.5989;
c_Bian = 0.0355;
d_Bian = -1.0456;
e_Bian = -9.5861*10^-4;

ln_y_Bian = zeros(4,6);
for n=1:4
    for i=1:6
        ln_y_Bian(n,i) = a_Bian + (b_Bian/T(n)) + ((c_Bian*rho_CO2(n,i))/T(n)) + (d_Bian + e_Bian*rho_CO2(n,i))*log(rho_CO2(n,i));
    end
end
% there is something wrong here
y_Bian = exp(ln_y_Bian);

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

figure(5);
% Subplot 1: Solubility vs. Pressure (S)
subplot(2,2,1);
hold on;
scatter(P, S(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r');  
scatter(P, S(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g');  
scatter(P, S(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b');  
scatter(P, S(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm');  
xlabel('Pressure [bar]');
ylabel('Drug solubility');
title('Solubility vs. Pressure (S)');
legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K');
xlim([115 280]);
ytickformat('%.2f');

% Subplot 2: Solubility vs. Pressure (y_b)
subplot(2,2,2);
hold on;
scatter(P, y_b(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r');
scatter(P, y_b(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g');  
scatter(P, y_b(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b');  
scatter(P, y_b(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm');  
xlabel('Pressure [bar]');
ylabel('Drug solubility');
title('Solubility vs. Pressure (y_b)');
legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K');
xlim([115 280]);
ytickformat('%.2f');
yticks(0.4:1:8.4);

% Subplot 3: Solubility vs. Density (S)
subplot(2,2,3);
hold on;
scatter(rho_CO2(1,:), S(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r');
scatter(rho_CO2(2,:), S(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g');
scatter(rho_CO2(3,:), S(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b');
scatter(rho_CO2(4,:), S(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm');
xlabel('Density [kg/m3]');
ylabel('Drug solubility');
title('Solubility vs. Density (S)');
legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K');
xlim([260 960]);

% Subplot 4: Solubility vs. Density (y_b)
subplot(2,2,4);
hold on;
scatter(rho_CO2(1,:), y_b(1,:), 20, 'r', 'o', 'MarkerFaceColor', 'r');
scatter(rho_CO2(2,:), y_b(2,:), 20, 'g', 's', 'MarkerFaceColor', 'g');
scatter(rho_CO2(3,:), y_b(3,:), 20, 'b', 'd', 'MarkerFaceColor', 'b');
scatter(rho_CO2(4,:), y_b(4,:), 20, 'm', '^', 'MarkerFaceColor', 'm');
xlabel('Density [kg/m3]');
ylabel('Drug solubility');
title('Solubility vs. Density (y_b)');
legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K');
xlim([260 960]);
ytickformat('%.2f');
yticks(0.4:1:8.4);
