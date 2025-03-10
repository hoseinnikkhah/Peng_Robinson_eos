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
rho_range = linspace(340,940,100);
% KJ Model info
a_KJ = -2.4761;
b_KJ = 5.2409*10^-3;
c_KJ = -4550.5306;

ln_y_KJ = zeros(4,6);
ln_y_KJ_cT = zeros(4,6);
for n=1:4
    for i=1:6
        ln_y_KJ(n,i) = a_KJ + b_KJ*rho_CO2(n,i) + c_KJ/T(n);
        ln_y_KJ_cT(n,i) = ln_y_KJ(n,i) - c_KJ/T(n);
    end
end

figure(1);
hold on;
plot(rho_CO2(1,:),ln_y_KJ_cT(1,:));
plot(rho_CO2(2,:),ln_y_KJ_cT(2,:));
plot(rho_CO2(3,:),ln_y_KJ_cT(3,:));
plot(rho_CO2(4,:),ln_y_KJ_cT(4,:));

legend()
xlabel('Density (kg/m^3)');
ylabel('lny - c/T');

figure(2);
hold on;
plot(rho_CO2(1,:),ln_y_KJ(1,:));
plot(rho_CO2(2,:),ln_y_KJ(2,:));
plot(rho_CO2(3,:),ln_y_KJ(3,:));
plot(rho_CO2(4,:),ln_y_KJ(4,:));

legend()
xlabel('Density (kg/m^3)');
ylabel('lny - c/T');
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

