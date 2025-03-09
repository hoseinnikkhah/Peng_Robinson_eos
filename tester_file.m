% PSub_1 calculations
PSub_1 = zeros(1,4);
deltaH_vap = 56.91;             % [kJ/mol]
DeltaH_vap = deltaH_vap * 10^3; % [J/mol]
T_ref = 304.18;                 % [K]
P_ref = 10^-6;                  % [bar]
for temp=1:4
    PSub_1(temp) = P_ref * exp((DeltaH_vap/R)*(1/T_ref - 1/T(temp)));
end
