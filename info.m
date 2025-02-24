T_c = 253;
Mw_drug = 598.5;
Mw_CO2 = 44.01;

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

figure(1);
hold on;
scatter(P, S(1,:), 'r', 'o', 'MarkerFaceColor', 'r');  
scatter(P, S(2,:), 'g', 's', 'MarkerFaceColor', 'g');  
scatter(P, S(3,:), 'b', 'd', 'MarkerFaceColor', 'b');  
scatter(P, S(4,:), 'm', '^', 'MarkerFaceColor', 'm');  

xlabel('Pressure [bar]');
ylabel('Drug solubility');
title('The solubility data of Ceftriaxone sodium drug in SC-CO2');

legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K');
xlim([115 280]);
ytickformat('%.2f');

figure(2);
hold on;
scatter(P, y_b(1,:), 'r', 'o', 'MarkerFaceColor', 'r');  
scatter(P, y_b(2,:), 'g', 's', 'MarkerFaceColor', 'g');  
scatter(P, y_b(3,:), 'b', 'd', 'MarkerFaceColor', 'b');  
scatter(P, y_b(4,:), 'm', '^', 'MarkerFaceColor', 'm');  

xlabel('Pressure [bar]');
ylabel('Drug solubility');
title('The solubility data of Ceftriaxone sodium drug in SC-CO2');

legend('T = 308 K', 'T = 318 K', 'T = 328 K', 'T = 338 K');
xlim([115 280]);
ytickformat('%.2f');
yticks(0.4:1:8.4);

