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

% KJ Model info
a_KJ = -2.4761;
b_KJ = 5.2409*10^3;
c_KJ = -4550.5306;

ln_y = a_KJ + b_KJ/T(i) c_KJ*rho_drug

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

