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

% KJ Model info
a = -2.4761;
b = 5.2409*10^3;
c = -4550.5306;

% GM Model info
a = 1.8309;
b = -5710.4987;
c = 3.8342;

% Chrastil Model info
a = 4.5889;
b = -12.7501;
c = -6759.0921;

% Bartle et al. Model info
a = 13.2886;
b = 9.0827;
c = -6845.6966;

% Sung-Shim Model info
a = 5.7399;
b = -606.0322;
c = -447.5124;
d = -36.8120;

% Bian et al model info
a = 4.2910;
b = -4806.5989;
c = 0.0355;
d = -1.0456;
e = -9.5861*10^4;

% Sodeifian model info
a = 
b = 
c = 
d = 
e = 
f = 
